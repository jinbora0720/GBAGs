README for Application data
================

# South Korea (KR)

-   PM2.5 in KR from AirKorea
    -   File: /data/KR\_Feb\_pm25.rds
-   Processing
    -   /Raw Data/AirKorea/202002\_BJ.csv
    -   Station Information stored in /Raw
        Data/AirKorea/station\_list\_final\_BJ.csv
    -   Available at
        <https://www.airkorea.or.kr/web/detailViewDown?pMENU_NO=125>
    -   Run /Raw Data/AirKorea/make\_KR\_pm25.R to obtain final .rds
        file

# California (CA), USA

-   PM2.5 and distance to the nearest fire in CA
    -   File: /data/CA\_final\_20211028.RDS
    -   Variables
        -   `date`: date (2020-08-01 \~ 2020-10-22)
        -   `easting`: UTM conversion of `long` by crs = 26911
        -   `northing`: UTM conversion of `lat` by crs = 26911
        -   `pm25_mean`: pm2.5 daily mean (removed values ≤ 0), measured
            in *μ*g/*m*<sup>3</sup>
        -   `from`: source (either PurpleAir or EPA)
        -   `logpm25_mean`: log transformed pm2.5
        -   `dist_to_fire`: Euclidean distance to the nearest fire on a
            given day
    -   Run /Raw Data/make\_CA\_pm25.R to obtain final .RDS files
-   Fires in CA
    -   File: /data/CA\_fire\_20211028.RDS
    -   Variables
        -   `date`: date (2020-08-01 \~ 2020-10-22)
        -   `easting`: UTM conversion of `long` by crs = 26911
        -   `northing`: UTM conversion of `lat` by crs = 26911
    -   Run /Raw Data/make\_CA\_pm25.R to obtain final .RDS files
-   Raw data from three sources (EPA and PurpleAir for PM2.5, NOAA for
    fires) are used.

## 1. EPA

-   PM2.5 in CA from EPA
    -   /Raw Data/CA\_EPA\_pm25\_2020.csv
    -   Available at
        <https://www.epa.gov/outdoor-air-quality-data/download-daily-data>

## 2. PurpleAir

-   PM2.5 in CA from PurpleAir
    -   /Raw Data/CAplus\_daily\_summary\_20210126.rds
    -   Station Information stored in /Raw
        Data/CAplus\_stations\_metadata\_20201105.rds
    -   How to Download Raw data from PurpleAir:
        <https://docs.google.com/document/d/15ijz94dXJ-YAZLi9iZ_RaBwrZ4KtYeCy08goGBwnbCU/edit>
        -   `/home/data/PurpleAirdata/sensordata`: Michele runs his
            script and uploads data in the VM Rstudio from all the
            stations over designated dates.
        -   `/home/data/PurpleAirdata/sensordata_daily_summary`: I
            create daily summary files.
        -   `/home/data/PurpleAirdata/CAplus_daily_summary.rds`: I save
            rds file from all the stations located in CA.  
    -   Followed EPA correction
        -   If PA<sub>cf<sub>1</sub></sub> ≤ 343 *μ*g/*m*<sup>3</sup>,
            PM2.5 = 0.52PA<sub>cf<sub>1</sub></sub> − 0.086RH + 5.75
        -   Otherwise, PM2.5 =
            0.46PA<sub>cf<sub>1</sub></sub> + 0.000393PA<sub>cf<sub>1</sub></sub><sup>2</sup> + 2.97
        -   PA<sub>cf<sub>1</sub></sub>=`pm25b_mean`: CF=1
            *μ*g/*m*<sup>3</sup> (removed values ≤ 0)
        -   RH = `RH`: relative humidity in %

## 3. Fire

-   File: /Raw Data/CA\_fire\_20210201.rds
-   Original txt files available at
    <https://www.ospo.noaa.gov/Products/land/hms.html>
-   Run /Raw Data/make\_CA\_fire.R to get the rds file.
