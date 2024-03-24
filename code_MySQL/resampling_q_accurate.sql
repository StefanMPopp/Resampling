USE resampling;
-- Find 5 best fitting resampling method (& threshold) for each ground truth parameter set
SELECT * FROM(
    WITH means AS ( -- calculating accuracy in % diff from ground truth
        SELECT
            paramSet
            , method
            , threshNr
            , ABS(1-AVG(angleMu)) AS accuracy
        FROM results_rel
        WHERE paramSet < 33 -- 33+ don't have accuracy
        GROUP BY paramSet, method, threshNr
    )
    SELECT
        m.paramSet
        , method
        , threshNr
        , m.accuracy
    FROM means m    
    JOIN ( -- Getting values of best
        SELECT *
        FROM (SELECT
                paramSet
                , accuracy
                , ROW_NUMBER() OVER(PARTITION BY paramSet ORDER BY accuracy) as acc_rank
            FROM means
        ) x
        WHERE acc_rank < 6 -- Change number to however many best you want
    ) b
    ON m.paramSet = b.paramSet and m.accuracy = b.accuracy
) top_res
JOIN params -- to know what the ground truth params are
ON top_res.paramSet = params.set_id
