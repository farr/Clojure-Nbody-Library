(ns nbody.test.hermite
  (:import java.util.Random)
  (:use clojure.test nbody.vec nbody.ics nbody.energy
        nbody.hermite))

(deftest single-step-order
  (let [bs (hot-spherical (Random.) 10)
        e1 (energy (advance bs 2e-8 1.0))
        e2 (energy (advance bs 1e-8 1.0))
        r (Math/abs (double (/ (- e1 -0.25)
                               (- e2 -0.25))))]
    (println "e1 = " e1 " e2 = " e2 " e = " (energy bs))
    (println r)))