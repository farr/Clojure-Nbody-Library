(ns nbody.test.gravity
  (:use clojure.test nbody.vec nbody.gravity nbody.test.predicates))

(deftest acc-test
  (let [m1 0.64626
        m2 0.196248
        r1 (double-array [0.680493 0.289422 0.497264])
        r2 (double-array [0.361698 0.826252 0.0804898])
        a (double-array 3 0.0)
        j (double-array 3 0.0)]
    (acc-and-jerk m1 r1 (double-array 3 0.0) m2 r2 (double-array 3 0.0) a j)
    (is (vector-close? a [-0.147896 0.249047 -0.193351] {:epsrel 1e-3}))))

(deftest jerk-test
  (let [m1 0.64626
        m2 0.196248
        r1 (double-array [0.680493 0.289422 0.497264])
        r2 (double-array [0.361698 0.826252 0.0804898])
        v1 (double-array [0.387609 0.765979 0.248525])
        v2 (double-array [0.0913329 0.122995 0.0394023])
        a (double-array 3 0.0)
        j (double-array 3 0.0)]
    (acc-and-jerk m1 r1 v1 m2 r2 v2 a j)
    (is (vector-close? j [-0.266232 -0.0814312 -0.26538] {:epsrel 1e-3}))))