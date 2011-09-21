(ns nbody.energy
  (:use nbody.vec))

(defn kinetic-energy
  "Returns the kinetic energy of a system of bodies."
  [bs]
  (reduce
   (fn [ke b]
     (let [v (doubles (:v b))
           m (double (:m b))
           ke (double ke)]
       (+ ke (* 0.5
                (* m (dot v v))))))
   0.0
   bs))

(defn- body-potential
  "Returns the potential energy of the given body with respect to the array of bodies bs."
  [b bs]
  (let [r (doubles (:r b))
        m (double (:m b))]
    (reduce
     (fn [pe b2]
       (if (not (= b b2))
         (let [r2 (doubles (:r b2))
               m2 (double (:m b2))
               pe (double pe)]
           (- pe (/ (* m m2)
                    (distance r r2))))
         pe))
     0.0
     bs)))

(defn potential-energy
  "Returns the potential energy of a system of bodies."
  [bs]
  (* 0.5 (reduce + (pmap (fn [b] (body-potential b bs)) bs))))

(defn energy
  "Returns the total energy of the bs system."
  [bs]
  (+ (kinetic-energy bs) (potential-energy bs)))