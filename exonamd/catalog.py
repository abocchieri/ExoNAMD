import os
import numpy as np
import pandas as pd
import requests
from datetime import datetime
from datetime import timedelta
from loguru import logger

from exonamd.utils import ROOT


@logger.catch
def download_nasa_confirmed_planets(
    min_sy_pnum=1,
    from_scratch=False,
):
    logger.info("Downloading NASA Exoplanet Archive confirmed planets")
    if from_scratch:
        df_old = None
        latest = datetime.strptime("1990-01-01", "%Y-%m-%d")
    else:
        df_old = pd.read_csv(os.path.join(ROOT, "data", "exo.csv"))
        latest = df_old["rowupdate"].max()
        latest = datetime.strptime(latest, "%Y-%m-%d")
        latest = latest - timedelta(days=1)
    latest = latest.strftime("%Y-%m-%d")

    logger.debug("Defining the SQL query to retrieve the required data")
    query = f"""
    SELECT 
        hostname, 
        pl_name, 
        default_flag,
        rowupdate,
        sy_pnum, 
        st_rad,
        st_mass,
        pl_orbper,
        pl_orbsmax, 
        pl_orbsmaxerr1, 
        pl_orbsmaxerr2, 
        pl_rade,
        pl_radeerr1,
        pl_radeerr2,
        pl_bmasse, 
        pl_bmasseerr1, 
        pl_bmasseerr2, 
        pl_orbeccen, 
        pl_orbeccenerr1, 
        pl_orbeccenerr2, 
        pl_orbincl, 
        pl_orbinclerr1, 
        pl_orbinclerr2,
        pl_trueobliq,
        pl_trueobliqerr1,
        pl_trueobliqerr2,
        pl_ratdor,
        pl_ratror
    FROM ps
    WHERE
        sy_pnum >= '{min_sy_pnum}'
        AND rowupdate > '{latest}'
    """

    logger.debug("Making the request to the API")
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
    params = {
        "query": query,
        "format": "json",
    }
    response = requests.get(url, params=params)

    if response.status_code != 200:
        logger.error(f"Error: {response.status_code} in fetching data")
        raise ValueError(f"Error: {response.status_code} in fetching data")

    data = response.json()
    df = pd.DataFrame(data)
    df = df.replace({None: np.nan, "": np.nan})

    logger.info("Data fetched")

    return df, df_old
