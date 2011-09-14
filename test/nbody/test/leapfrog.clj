(ns nbody.test.leapfrog
  (:import java.util.Random)
  (:use clojure.test nbody.leapfrog nbody.ics
        nbody.energy))

(deftest leapfrog-energy-scaling
  (let [bs (map map->particle (hot-spherical (Random.) 10))
        e1 (energy (advance bs 1e-5))
        e2 (energy (advance bs 1e-6))
        r (Math/abs (double (/ (- e1 -0.25)
                               (- e2 -0.25))))]
    (is (and (< r 11.31) ; Geometric mean of 8 and 16
             (> r 5.66))))) ; Geometric mean of 4 and 8