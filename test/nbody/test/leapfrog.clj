(ns nbody.test.leapfrog
  (:import java.util.Random)
  (:use clojure.test nbody.leapfrog nbody.ics
        nbody.energy nbody.vec))

(deftest leapfrog-energy-scaling
  (let [bs (map map->particle (hot-spherical (Random.) 100))
        e1 (energy (advance bs 1e-3))
        e2 (energy (advance bs 5e-4))
        r (Math/abs (double (/ (- e1 -0.25)
                               (- e2 -0.25))))]
    (is (and (< r 11.31) ; Geometric mean of 8 and 16
             (> r 5.66))))) ; Geometric mean of 4 and 8

(deftest long-term-leapfrog-scaling
  (let [bs (map map->particle (hot-spherical (Random.) 100))]
    (let [bs1 (nth (iterate #(advance % 1e-5) bs) 100)
          bs2 (nth (iterate #(advance % 5e-6) bs) 200)]
      (let [e1 (energy bs1)
            e2 (energy bs2)]
        (let [r (/ (+ e1 0.25)
                   (+ e2 0.25))]
          (is (and (< r (Math/sqrt (* 4.0 8.0)))
                   (> r (Math/sqrt (* 2.0 4.0))))))))))

(deftest circular-binary-test
  (let [bs (map map->particle
                (list {:id 1 :m 2 :t 0 :r (double-array [1.0 0.0 0.0]) :v (double-array [0.0 0.0 (Math/sqrt 0.5)])}
                      {:id 2 :m 2 :t 0 :r (double-array [-1.0 0.0 0.0]) :v (double-array [0.0 0.0 (- (Math/sqrt 0.5))])}))
        new-bs1 (advance bs 0.3)
        new-bs2 (advance bs 0.15)]
    (let [e (energy bs)
          e1 (energy new-bs1)
          e2 (energy new-bs2)]
      (let [r (Math/abs (double (/ (- e e1) (- e e2))))]
        (is (> r (Math/sqrt (* 4.0 8.0))))))))

(deftest long-circular-binary-test
  (let [bs (map map->particle
                (list {:id 1 :m 1 :t 0 :r (double-array [0.0 0.0 0.5]) :v (double-array [0.0 (Math/sqrt 0.5) 0.0])}
                      {:id 2 :m 1 :t 0 :r (double-array [0.0 0.0 -0.5]) :v (double-array [0.0 (- (Math/sqrt 0.5)) 0.0])}))]
    (let [new-bs1 (nth (iterate #(advance % 0.01) bs) 100)
          new-bs2 (nth (iterate #(advance % 0.005) bs) 200)]
       (let [e (energy bs)
             e1 (energy new-bs1)
             e2 (energy new-bs2)]
         (let [r (Math/abs (double (/ (- e e1) (- e e2))))]
           (is (and (> r (Math/sqrt (* 2.0 4.0)))))
           (is (< (Math/abs (double (- 1.0 (distance (:r (first new-bs1))
                                                     (:r (first (rest new-bs1)))))))
                  1e-3)))))))

(deftest leapfrog-eight
  (let [r (double-array [-0.97000436 0.24308753 0.0])
        v (double-array [0.93240737 0.86473146 0.0])]
    (let [bs (map map->particle
                  (list {:id 1 :m 1 :t 0 :r r :v (v* -0.5 v)}
                        {:id 2 :m 1 :t 0 :r (v* -1.0 r) :v (v* -0.5 v)}
                        {:id 3 :m 1 :t 0 :r (double-array 3 0.0) :v v}))]
      (let [new-bs1 (nth (iterate #(advance % 0.001) bs) 1000)
            new-bs2 (nth (iterate #(advance % 0.0005) bs) 2000)]
        (let [e (energy bs)
              e1 (energy new-bs1)
              e2 (energy new-bs2)]
          (let [r (Math/abs (double (/ (- e e1) (- e e2))))]
            (is (and (< r (Math/sqrt (* 4.0 8.0)))
                     (> r (Math/sqrt (* 2.0 4.0)))))))))))