(ns nbody.test.ics
  (:import java.util.Random)
  (:use clojure.test nbody.vec nbody.ics nbody.energy nbody.test.predicates))

(deftest hot-spherical-virial
  (let [bs (hot-spherical (Random.) 100)
        Q (/ (kinetic-energy bs)
             (- (potential-energy bs)))]
    (is (and (> Q 0.4)
             (< Q 0.6)))))

(deftest hot-spherical-energy
  (let [bs (hot-spherical (Random.) 100)]
    (is (close? (energy bs) -0.25))))

(deftest hot-spherical-com
  (let [com (com (hot-spherical (Random.) 100))]
    (is (close? (norm com) 0.0))))

(deftest hot-spherical-vtot
  (let [vtot (vtot (hot-spherical (Random.) 100))]
    (is (close? (norm vtot) 0.0))))