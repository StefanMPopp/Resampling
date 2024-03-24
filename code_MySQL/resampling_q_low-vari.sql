USE resampling;
-- Find 5 best fitting resampling method (& threshold) for each ground truth parameter set
-- and select of those the one with the lowest variablity
WITH acc_var AS(
    SELECT * FROM(
        WITH acc AS(
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
                ON m.paramSet = b.paramSet AND m.accuracy = b.accuracy
            ) top_res
            JOIN params -- to know what the ground truth params are
            ON top_res.paramSet = params.set_id
        )
        SELECT
            a.paramSet
            , a.method
            , a.threshNr
            , a.accuracy
            , a.length
            , a.angle_mu
            , v.rnge95_angleMu
        FROM acc a
        JOIN results_var v
        ON a.paramSet = v.paramSet AND a.method = v.method AND a.threshNr = v.threshNr
    ) joined
)
SELECT
    av.paramSet
    , av.method
    , av.threshNr
    , av.accuracy
    , av.rnge95_angleMu
    , min_var
FROM acc_var av
JOIN
    (
    SELECT
        paramSet
        , MIN(rnge95_angleMu) AS min_var
    FROM acc_var
    GROUP BY paramSet
    ) mn
ON av.paramSet = mn.paramSet AND av.rnge95_angleMu = mn.min_var