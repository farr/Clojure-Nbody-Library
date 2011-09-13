(ns nbody.test.predicates
  (use clojure.test))

(defn close?
  "True only if the difference between x and y satisfies the given
  absolute or relative error bounds."
  ([x y]
     (close? x y {:epsabs 1e-8 :epsrel 1e-8}))
  ([x y {:keys [epsabs epsrel]}]
     (let [eabs (double (or epsabs 1e-8))
           erel (double (or epsrel 1e-8))
           x (double x)
           y (double y)]
       (let [dx (Math/abs (- x y))
             ave (* 0.5 (+ (Math/abs x) (Math/abs y)))]
         (<= dx (+ eabs (* ave erel)))))))

(defn vector-close?
  "True only if each of the components of x and y are close?."
  ([x y]
     (vector-close? x y {:epsabs 1e-8 :epsrel 1e-8}))
  ([x y {:keys [epsabs epsrel] :as e}]
     (loop [x (seq x) y (seq y)]
       (or (not (seq x))
           (not (seq y))
           (and (close? (first x) (first y) e)
                (recur (rest x) (rest y)))))))