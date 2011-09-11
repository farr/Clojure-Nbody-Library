(ns nbody.test.vec
  (:use clojure.test nbody.vec nbody.test.predicates))

(deftest dot-test
  (let [x (double-array [1 2 3])
        y (double-array [2 3 4])]
    (is (close? (dot x y) 20.0))))

(deftest distance-test
  (let [x (double-array [1 2 3])
        y (double-array [2 3 4])]
    (is (close? (distance x y) (Math/sqrt 3)))))