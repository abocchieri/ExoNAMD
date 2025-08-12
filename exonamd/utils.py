import os
import urllib
import concurrent.futures
import re
import pickle
import numpy as np
import pandas as pd
import requests
import warnings
from astroquery.simbad import Simbad as simbad
from loguru import logger


ROOT = os.path.realpath(os.path.dirname(__file__))


# Functions below are originally from gen_tso
# by Patricio Cubillos: https://github.com/pcubillos/gen_tso
# modified to work in this project


def get_host(name):
    """
    Extract host name from a given planet name.
    Valid confirmed planet names end with a lower-case letter preceded
    by a blank.  Valid planet candidate names end with a dot followed
    by two numbers.

    Examples
    --------
    >>> get_host('TOI-741.01')
    >>> get_host('WASP-69 b')
    """
    if is_letter(name):
        return name[:-2]
    if "." in name:
        idx = name.rfind(".")
        return name[:idx]
    return ""


def get_letter(name):
    """
    Extract 'letter' identifier for a planet name.
    Valid confirmed planet names end with a lower-case letter preceded
    by a blank.  Valid planet candidate names end with a dot followed
    by two numbers.

    Examples
    --------
    >>> get_letter('TOI-741.01')
    >>> get_letter('WASP-69 b')
    """
    if is_letter(name):
        return name[-2:]
    if "." in name:
        idx = name.rfind(".")
        return name[idx:]
    return ""


def is_letter(name):
    """
    Check if name ends with a blank + lower-case letter (it's a planet)
    """
    return name[-1].islower() and name[-2] == " "


def is_candidate(name):
    """
    Check if name ends with a blank + lower-case letter (it's a planet)
    """
    return len(name) >= 3 and name[-3] == "." and name[-2:].isnumeric()


def invert_aliases(aliases):
    """
    Invert an {alias:name} dictionary into {name:aliases_list}
    """
    aka = {}
    for key, val in aliases.items():
        if val not in aka:
            aka[val] = []
        aka[val].append(key)
    return aka


def fetch_nea_aliases(targets):
    """
    Fetch target aliases as known by https://exoplanetarchive.ipac.caltech.edu/

    Note 1: a search of a planet or stellar target returns the
        aliases for all bodies in that planetary system.
    Note 2: there might be more than one star per system

    Parameters
    ----------
    targets: String or 1D iterable of strings
        Target(s) to fetch from the NEA database.

    Returns
    -------
    host_aliases_list: 1D list of dictionaries
        List of host-star aliases for each target.
    planet_aliases_list: 1D list of dictionaries
        List of planetary aliases for each target.

    Examples
    --------
    >>> import gen_tso.catalogs as cat
    >>> targets = ['WASP-8 b', 'KELT-7', 'HD 189733']
    >>> host_aliases, planet_aliases = cat.fetch_nea_aliases(targets)

    >>> host_aliases, planet_aliases = cat.fetch_nea_aliases('WASP-999')
    """
    if isinstance(targets, str):
        targets = [targets]
    ntargets = len(targets)

    urls = np.array(
        [
            "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/Lookup/"
            f"nph-aliaslookup.py?objname={urllib.parse.quote(target)}"
            for target in targets
        ]
    )

    def fetch_url(url):
        try:
            response = requests.get(url)
            return response
        except:
            return None

    fetch_status = np.tile(2, ntargets)
    responses = np.tile({}, ntargets)
    n_attempts = 0
    while np.any(fetch_status > 0) and n_attempts < 10:
        n_attempts += 1
        mask = fetch_status > 0
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(fetch_url, urls[mask]))

        j = 0
        for i in range(ntargets):
            if fetch_status[i] <= 0:
                continue
            r = results[j]
            j += 1
            if r is None:
                continue
            if not r.ok:
                warnings.warn(f"Alias fetching failed for '{targets[i]}'")
                fetch_status[i] -= 1
                continue
            responses[i] = r.json()
            fetch_status[i] = 0
        fetched = np.sum(fetch_status <= 0)
        logger.info(f"Fetched {fetched}/{ntargets} entries on try {n_attempts}")

    host_aliases_list = []
    planet_aliases_list = []
    for i, resp in enumerate(responses):
        if resp == {}:
            logger.warning(f"NEA alias fetching failed for '{targets[i]}'")
            host_aliases_list.append({})
            planet_aliases_list.append({})
            continue
        if resp["manifest"]["lookup_status"] == "System Not Found":
            logger.warning(f"NEA alias not found for '{targets[i]}'")
            host_aliases_list.append({})
            planet_aliases_list.append({})
            continue

        host_aliases = {}
        star_set = resp["system"]["objects"]["stellar_set"]["stars"]
        for star in star_set.keys():
            if "is_host" not in star_set[star]:
                continue
            for alias in star_set[star]["alias_set"]["aliases"]:
                host_aliases[alias] = star
        host_aliases_list.append(host_aliases)

        planet_aliases = {}
        planet_set = resp["system"]["objects"]["planet_set"]["planets"]
        for planet in planet_set.keys():
            for alias in planet_set[planet]["alias_set"]["aliases"]:
                planet_aliases[alias] = planet
        planet_aliases_list.append(planet_aliases)

    return host_aliases_list, planet_aliases_list


def get_children(host_aliases, planet_aliases):
    """
    Cross check a dictionary of star and planet aliases to see
    whether the star is the host of the planets.
    """
    # get all planet aliases minus the 'letter' identifier
    planet_aka = invert_aliases(planet_aliases)
    for planet, aliases in planet_aka.items():
        aka = []
        for alias in aliases:
            len_letter = len(get_letter(alias))
            aka.append(alias[0:-len_letter])
        planet_aka[planet] = aka

    # cross_check with host aliases
    children = []
    for planet, aliases in planet_aka.items():
        if np.any(np.in1d(aliases, host_aliases)):
            children.append(planet)

    aliases = {
        alias: planet for alias, planet in planet_aliases.items() if planet in children
    }
    return aliases


def fetch_simbad_aliases(target, verbose=True):
    """
    Fetch target aliases and Ks magnitude as known by Simbad.

    Examples
    --------
    >>> from gen_tso.catalogs.update_catalogs import fetch_simbad_aliases
    >>> aliases, ks_mag = fetch_simbad_aliases('WASP-69b')
    """
    simbad.reset_votable_fields()
    simbad.remove_votable_fields("coordinates")
    simbad.add_votable_fields("otype", "otypes", "ids")
    simbad.add_votable_fields("flux(K)")
    # simbad.add_votable_fields("fe_h")

    host_alias = []
    kmag = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        simbad_info = simbad.query_object(target)
    if simbad_info is None:
        if verbose:
            logger.warning(f"no Simbad entry for target {repr(target)}")
        return host_alias, kmag

    object_type = simbad_info["OTYPE"].value.data[0]
    if "Planet" in object_type:
        if target[-1].isalpha():
            host = target[:-1]
        elif "." in target:
            end = target.rindex(".")
            host = target[:end]
        else:
            # target_id = simbad_info['MAIN_ID'].value.data[0]
            return host_alias, kmag
        # go after star
        simbad_info = simbad.query_object(host)
        if simbad_info is None:
            if verbose:
                logger.warning(f"Simbad host {repr(host)} not found")
            return host_alias, kmag

    host_info = simbad_info["IDS"].value.data[0]
    host_alias = host_info.split("|")
    kmag = simbad_info["FLUX_K"].value.data[0]
    # fetch metallicity?
    if not np.isfinite(kmag):
        kmag = np.nan
    return host_alias, kmag


@logger.catch
def fetch_aliases(hosts, output_file=None, known_aliases=None):
    """
    Fetch aliases from the NEA and Simbad databases for a list
    of host stars.  Store output dictionary of aliases to pickle file.

    Parameters
    ----------
    hosts: List of strings
        Host star names of targets to search.
    output_file: String
        If not None, save outputs to file (as pickle binary format).
    known_aliases: Dictionary
        Dictionary of known aliases, the new aliases will be added
        on top of this dictionary.

    Returns
    -------
    aliases: Dictionary
        Dictionary of aliases with one entry per system where
        the key is the host name.  Each entry is a dictionary
        containing:
            'host' the host name (string)
            'planets': array of planets in the system
            'host_aliases': list of host aliases
            'planet_aliases': dictionary of planets' aliases with
                (key,value) as (alias_name, name)

    Examples
    --------
    >>> import gen_tso.catalogs as cat
    >>> from gen_tso.utils import ROOT

    >>> # Confirmed targets
    >>> targets = cat.load_targets()
    >>> hosts = np.unique([target.host for target in targets])
    >>> output_file = f'{ROOT}data/nea_aliases.pickle'
    >>> aliases = cat.fetch_aliases(hosts, output_file)
    """

    logger.info("Fetching aliases for host stars and their planets")

    if known_aliases is None:
        known_aliases = {}

    host_aliases, planet_aliases = fetch_nea_aliases(hosts)

    aliases = {}
    nhosts = len(hosts)
    for i in range(nhosts):
        # Isolate host-planet(s) aliases
        stars = np.unique(list(host_aliases[i].values()))
        hosts_aka = invert_aliases(host_aliases[i])
        for host, h_aliases in hosts_aka.items():
            if hosts[i] in h_aliases:
                host_name = host
                break

        h_aliases = [
            alias for alias, host in host_aliases[i].items() if host == host_name
        ]
        if len(stars) == 1:
            p_aliases = planet_aliases[i].copy()
        else:
            p_aliases = get_children(h_aliases, planet_aliases[i])
        p_aliases = {
            re.sub(r"\s+", " ", key): val
            for key, val in p_aliases.items()
            if is_letter(key) or is_candidate(key)
        }
        children_names = np.unique(list(p_aliases.values()))

        # Complement with Simbad aliases:
        s_aliases, kmag = fetch_simbad_aliases(host_name, verbose=False)
        new_aliases = []
        for alias in s_aliases:
            alias = re.sub(r"\s+", " ", alias)
            is_new = (
                alias.startswith("G ")
                or alias.startswith("GJ ")
                or alias.startswith("Wolf ")
                or alias.startswith("2MASS ")
            )
            if is_new and alias not in h_aliases:
                new_aliases.append(alias)
                h_aliases.append(alias)

        # Replicate host aliases as planet aliases:
        planet_aka = invert_aliases(p_aliases)
        for planet, pals in planet_aka.items():
            for host in h_aliases:
                letter = get_letter(planet)
                planet_name = f"{host}{letter}"
                # The conditions to add a target:
                is_new = planet_name not in pals
                # There is a planet or a candidate in list
                planet_exists = np.any(
                    [get_host(p) == host and is_letter(p) for p in pals]
                )
                candidate_exists = np.any(
                    [get_host(p) == host and is_candidate(p) for p in pals]
                )
                # Do not downgrade planet -> candidate
                not_downgrade = not (is_candidate(planet_name) and planet_exists)
                # No previous alias (hold-off TESS names)
                new_entry = (
                    not planet_exists
                    and not candidate_exists
                    and not planet_name.startswith("TOI")
                )
                # There is a letter version of it with same root
                letter_exists = np.any(
                    [p.startswith(host) and is_letter(p) for p in pals]
                )
                # Upgrade candidate->planet only if is lettered anywhere else
                upgrade = is_letter(planet_name) and candidate_exists and letter_exists
                if is_new and not_downgrade and (new_entry or upgrade):
                    p_aliases[planet_name] = planet

        system = {
            "host": host_name,
            "planets": children_names,
            "host_aliases": h_aliases,
            "planet_aliases": p_aliases,
        }
        aliases[host_name] = system

    # Add previously known aliases (but give priority to the new ones)
    for host in list(known_aliases):
        if host not in aliases:
            aliases[host] = known_aliases[host]

    logger.info("Aliases fetched.")

    if output_file is not None:
        logger.info(f"Saving aliases to pickle")
        with open(output_file, "wb") as handle:
            pickle.dump(aliases, handle, protocol=4)
        logger.info(f"Aliases saved to {output_file}")

    return aliases


@logger.catch
def update_host(row, aliases):
    host = get_value(row["hostname"])
    for key, item in aliases.items():
        if host in item["host_aliases"]:
            if host != key:
                logger.debug(f"Updating host {host} to {key}")
            return key
    return host


@logger.catch
def update_planet(row, aliases):
    planet = get_value(row["pl_name"])
    for key, item in aliases.items():
        planet_aliases = item["planet_aliases"]
        if planet in planet_aliases.keys():
            name = planet_aliases[planet]
            if name[: len(key)] != key:
                logger.debug(f"Updating planet {planet} to {key + name[-2:]}")
                return key + name[-2:]
            elif planet != name:
                logger.debug(f"Updating planet {planet} to {name}")
                return name
    return planet


@logger.catch
def groupby_apply_merge(df, groupby, func, *args, allow_overwrite=False, **kwargs):
    retval = df.groupby(groupby).apply(func, *args, **kwargs)

    if allow_overwrite:
        # Identify columns that would cause a conflict during the merge
        overlapping_columns = retval.columns.intersection(df.columns)

        # Drop those columns from the original dataframe
        df = df.drop(columns=overlapping_columns)

    # Perform the merge
    return df.merge(retval, left_on=groupby, right_index=True)


@logger.catch
def check_name(names):
    if len(set(name[:3] for name in names)) > 1:
        return False
    return True


def get_value(value):
    return value.iloc[0] if isinstance(value, pd.Series) else value


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import truncnorm, norm


def sample_trunc_normal(
    mu=0.0, sigma=1.0, lower=-1.0, upper=2.0, n=10000, random_state=None
):
    """
    Sample from a truncated normal distribution using scipy.stats.truncnorm.

    Parameters
    ----------
    mu, sigma : float
        Mean and standard deviation of the underlying normal.
    lower, upper : float
        Truncation bounds (inclusive).
    n : int
        Number of samples to draw.
    random_state : int | None
        Seed for reproducibility.

    Example
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.stats import truncnorm, norm

    >>> mu = 1.0
    >>> sigma = 0.8
    >>> lower = 0.0
    >>> upper = 2.0
    >>> n = 10000
    >>> seed = 42

    >>> samples = sample_trunc_normal(mu=mu, sigma=sigma, lower=lower, upper=upper, n=n, random_state=seed)

    >>> # Plot histogram and theoretical PDF
    >>> fig, ax = plt.subplots(figsize=(8, 4.5))
    >>> ax.hist(samples, bins=60, density=True, alpha=0.6)
    >>> ax.set_xlabel("x")
    >>> ax.set_ylabel("Density")
    >>> ax.set_title(f"Truncated Normal: mu={mu}, sigma={sigma}, lower={lower}, upper={upper}")

    >>> # Theoretical PDF
    >>> Z = norm.cdf((upper - mu) / sigma) - norm.cdf((lower - mu) / sigma)
    >>> xs = np.linspace(lower, upper, 400)
    >>> pdf = norm.pdf((xs - mu) / sigma) / (sigma * Z)
    >>> ax.plot(xs, pdf, linewidth=2)

    >>> plt.show()

    Returns
    -------
    samples : ndarray, shape (n,)
    """
    if random_state is not None:
        np.random.seed(random_state)

    # Convert bounds to standard normal units
    if sigma == 0:
        sigma = 1e-9
    a, b = (lower - mu) / sigma, (upper - mu) / sigma
    return truncnorm.rvs(a, b, loc=mu, scale=sigma, size=n)
