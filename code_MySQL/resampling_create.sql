-- CREATE DATABASE resampling
USE resampling;

CREATE Table params(
    length INT,
    step_distr DECIMAL,
    step_mu DECIMAL,
    step_sd DECIMAL,
    angle_mu DECIMAL,
    angle_sd DECIMAL,
    noise_angle DECIMAL,
    noise_xy DECIMAL,
    smooth_w INT
);

-- Import all .csv
-- Parameters for making ground truth tracks
SET GLOBAL local_infile=1;
LOAD DATA LOCAL INFILE '/home/stefan/Documents/SQL/resampling/params.csv'
INTO TABLE params
FIELDS TERMINATED BY ',' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;

ALTER TABLE params  -- add PRIMARY
ADD COLUMN set_id INT AUTO_INCREMENT
PRIMARY KEY FIRST;

-- -- Methods (not useful?)
-- CREATE TABLE methods(
--     name VARCHAR(20) PRIMARY
-- )

-- Thresholds
CREATE TABLE thresholds(
    method VARCHAR(20),
    thresh1 DECIMAL,
    thresh2 DECIMAL,
    thresh DECIMAL,
    threshNr INT
);
LOAD DATA LOCAL INFILE '/home/stefan/Documents/SQL/resampling/threshs.csv'
INTO TABLE thresholds
FIELDS TERMINATED BY ',' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;

ALTER TABLE thresholds
DROP COLUMN thresh1,
DROP thresh2;

-- Resampled data absolute
CREATE Table results_abs(
    paramSet INT,
    rep INT,
    method varchar(20),
    threshNr INT,
    nrTurns INT,
    length INT,
    angleMu DECIMAL(10,2),
    angleSD DECIMAL(10,2),
    stepMu DECIMAL(10,2),
    stepSD DECIMAL(10,4)
);
LOAD DATA LOCAL INFILE '/home/stefan/Documents/SQL/resampling/results_abs.csv'
INTO TABLE results_abs
FIELDS TERMINATED BY ',' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;

-- Resampled data relative to ground truth
CREATE Table results_rel(
    res_id INT NOT NULL PRIMARY KEY,
    paramSet INT,
    rep INT,
    method varchar(20),
    threshNr INT,
    nrTurns INT,
    length DECIMAL(10,4),
    angleMu DECIMAL(10,4),
    angleSD DECIMAL(10,4),
    stepMu DECIMAL(10,4),
    stepSD DECIMAL(10,4)
);
LOAD DATA LOCAL INFILE '/home/stefan/Documents/SQL/resampling/results_rel.csv'
INTO TABLE results_rel
FIELDS TERMINATED BY ',' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;

-- resampled data variances between replications
CREATE Table results_var(
    paramSet INT,
    method VARCHAR(20),
    threshNr INT,
    rnge95_nrTurns DECIMAL(10,4),
    rnge95_length DECIMAL(10,4),
    rnge95_angleMu DECIMAL(10,4),
    rnge95_angleSD DECIMAL(10,4),
    rnge95_stepMu DECIMAL(10,4),
    rnge95_stepSD DECIMAL(10,4)
);
LOAD DATA LOCAL INFILE '/home/stefan/Documents/SQL/resampling/results_var.csv'
INTO TABLE results_var
FIELDS TERMINATED BY ',' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;