(ns nbody.body
  (:use nbody.vec))

(defrecord Body [id ^double t ^double m ^doubles x ^doubles v])

(defn ke
  "Returns the kinetic energy of the given body."
  [^Body b]
  (* (* 0.5 (.m b))
     (norm (.v b))))

(defn pe
  "Returns the potential energy between b1 and b2."
  [^Body b1 ^Body b2]
  (/ (* (.m b1) (.m b2))
     (distance (.x b1) (.x b2))))

(defn drift
  "Drifts the given body forward by dt."
  [^Body b dt]
  (let [dt (double dt)
        x (doubles (.x b))
        v (doubles (.v b))]
    (assoc b :x (amap x i xnew (+ (aget x i) (* (aget v i) dt))))))