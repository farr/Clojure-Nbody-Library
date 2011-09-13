(ns nbody.test.hermite
  (:import java.util.Random)
  (:use clojure.test nbody.vec nbody.ics nbody.energy
        nbody.hermite))

(deftest single-step-order
  (let [bs (hot-spherical (Random.) 100)
        e1 (energy (advance bs 1e-6 1.0))
        e2 (energy (advance bs 1e-7 1.0))
        r (Math/abs (double (/ (- e1 -0.25)
                               (- e2 -0.25))))]
    (println r)))