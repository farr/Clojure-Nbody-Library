(ns nbody.leapfrog
  (:use nbody.vec nbody.gravity))

(defrecord Particle [id ^double m ^double t
                     ^doubles r ^doubles v])

(defn map->particle
  "Turns a map into a leapfrog particle."
  [m]
  (if (instance? Particle m)
    m
    (new Particle
         (:id m)
         (:m m)
         (:t m)
         (:r m)
         (:v m))))

(defn drift
  "Returns a new set of particles, ps, drifted forward by dt."
  [ps dt]
  (pmap
   (fn [^Particle p]
     (let [r (doubles (.r p))
           v (doubles (.v p))
           t (double (.t p))
           dt (double dt)]
       (let [rnew (amap r i rnew (+ (aget r i) (* dt (aget v i))))]
         (assoc p :r rnew :t (+ t dt)))))
   ps))

(defn kick
  "Returns a new set of particles, ps, kicked by dt."
  [ps dt]
  (pmap
   (fn [^Particle p]
     (let [r (doubles (.r p))
           v (doubles (.v p))
           vnew (doubles (aclone v))
           a (doubles (double-array 3 0.0))
           dt (double dt)]
       (doseq [^Particle p2 ps]
         (when (not (= p p2))
           (let [r2 (doubles (.r p2))
                 m2 (double (.m p2))]
             (acc! r m2 r2 a)
             (dotimes [i 3]
               (aset vnew i (+ (aget vnew i) (* dt (aget a i))))))))
       (assoc p :v vnew)))
   ps))

(defn advance
  "Advances a system of bodies, bs, forward by dt."
  [bs dt]
  (drift (kick (drift bs (/ dt 2.0)) dt) (/ dt 2.0)))