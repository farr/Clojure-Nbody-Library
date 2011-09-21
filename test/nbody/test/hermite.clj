(ns nbody.test.hermite
  (:import (java.util Random))
  (:use clojure.test nbody.vec nbody.ics nbody.energy
        nbody.hermite nbody.test.predicates))

(deftest predict!-test
  (let [p (new nbody.hermite.Particle
               "id" 0.0 0.0 0.0
               (double-array [0.220107 0.887105 0.362987])
               (double-array [0.231283 0.388615 0.948466])
               (double-array [0.198629 0.602694 0.771291])
               (double-array [0.779636 0.181176 0.787191]))
        tp 0.71329
        rp (double-array 3 0.0)
        vp (double-array 3 0.0)]
    (predict! p tp rp vp)
    (is (vector-close? rp (double-array [0.482765 1.32858 1.28334]) {:epsrel 1e-3}))
    (is (vector-close? vp (double-array [0.571296 0.8646 1.69887]) {:epsrel 1e-3}))))

(deftest total-acc-and-jerk-test
  (let [ps (map map->particle (hot-spherical (Random.) 10))
        [a j] (reduce (fn [[a j] [aa jj]] [(v+ aa a) (v+ j jj)])
                      (map #(total-acc-and-jerk % ps) ps))]
    (is (vector-close? a (double-array 3 0.0)))
    (is (vector-close? j (double-array 3 0.0)))))

(deftest circular-binary-test
  (let [bs (setup
            (list {:id 1 :m 1 :t 0 :r (double-array [0.5 0.0 0.0]) :v (double-array [0.0 0.0 (Math/sqrt 0.5)])}
                  {:id 2 :m 1 :t 0 :r (double-array [-0.5 0.0 0.0]) :v (double-array [0.0 0.0 (- (Math/sqrt 0.5))])})
            1e-3)
        new-bs1 (map #(advanced-particle (assoc % :tnext 0.3) bs 1e-3) bs)
        new-bs2 (map #(advanced-particle (assoc % :tnext 0.15) bs 1e-3) bs)]
    (let [e (energy bs)
          e1 (energy new-bs1)
          e2 (energy new-bs2)]
      (let [r (Math/abs (double (/ (- e e1) (- e e2))))]
        (is (> r (Math/sqrt (* 16.0 32.0))))))))

(deftest long-circular-binary-test
  (let [bs (setup
            (list {:id 1 :m 1 :t 0 :r (double-array [0.0 0.0 0.5]) :v (double-array [0.0 (Math/sqrt 0.5) 0.0])}
                  {:id 2 :m 1 :t 0 :r (double-array [0.0 0.0 -0.5]) :v (double-array [0.0 (- (Math/sqrt 0.5)) 0.0])})
            1e-3)]
    (letfn [(advancer [dt]
              (fn [bs]
                (map #(advanced-particle (assoc % :tnext (+ (:t %) dt)) bs 1e-3) bs)))]
      (let [new-bs1 (nth (iterate (advancer 0.1) bs) 100)
            new-bs2 (nth (iterate (advancer 0.05) bs) 200)]
       (let [e (energy bs)
             e1 (energy new-bs1)
             e2 (energy new-bs2)]
         (let [r (Math/abs (double (/ (- e e1) (- e e2))))]
           (is (> r (Math/sqrt (* 8.0 16.0))))
           (is (< (Math/abs (double (- 1.0 (distance (:r (first new-bs1))
                                                     (:r (first (rest new-bs1)))))))
                  1e-3))))))))

(deftest advanced-body-test
  (let [ps (setup (hot-spherical (Random.) 10) 1e-3)]
    (let [new-ps1 (map #(advanced-particle (assoc % :tnext 0.01) ps 1e-3) ps)
          new-ps2 (map #(advanced-particle (assoc % :tnext 0.005) ps 1e-3) ps)]
      (let [e (energy ps)
            e1 (energy new-ps1)
            e2 (energy new-ps2)]
        (let [r (Math/abs (double (/ (- e1 e) (- e2 e))))]
          (is (and (< r (Math/sqrt (* 32.0 64.0)))
                   (> r (Math/sqrt (* 16.0 32.0))))))))))

(deftest finish-test
  (let [bs (setup (hot-spherical (Random.) 100) 1e-3)
        e (energy bs)]
    (let [e1 (energy (finish bs 1e-5 1e-3))
          e2 (energy (finish bs 5e-4 1e-3))]
     (let [err1 (Math/abs (double (/ (- e1 e) e)))
           err2 (Math/abs (double (/ (- e2 e) e)))]
       (let [r (/ err1 err2)]
         true)))))