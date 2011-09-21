(ns nbody.ics
  (import java.util.Random)
  (:use nbody.vec nbody.energy))

(defn mtot
  "Returns the total mass in a system of bodies."
  [bs]
  (reduce (fn [mtot b] (+ mtot (:m b))) 0.0 bs))

(defn com
  "Returns the center of mass of the bodies."
  [bs]
  (let [mtot (double (mtot bs))]
    (reduce
     (fn [com b]
       (let [r (doubles (:r b))
             m (double (:m b))
             com (doubles com)]
         (amap com i new-com
               (+ (aget com i)
                  (* (aget r i)
                     (/ m mtot))))))
     (double-array 3 0.0)
     bs)))

(defn vtot
  "Returns the center-of-mass velocity of the array of bodies."
  [bs]
  (let [mtot (double (mtot bs))]
    (reduce
     (fn [vtot b]
       (let [v (doubles (:v b))
             m (double (:m b))
             vtot (doubles vtot)]
         (amap vtot i new-vtot
               (+ (aget vtot i)
                  (* (aget v i)
                     (/ m mtot))))))
     (double-array 3 0.0)
     bs)))

(defn to-com-frame
  "Returns the system bs adjusted to the center-of-mass frame."
  [bs]
  (let [com (doubles (com bs))
        vtot (doubles (vtot bs))]
    (map (fn [b]
           (let [r (doubles (:r b))
                 v (doubles (:v b))]
             (let [rnew (v- r com)
                   vnew (v- v vtot)]
               (assoc b :r rnew :v vnew))))
         bs)))

(defn to-standard-units
  "Returns a new system of bodies which is in the standard units of
  Heggie and Mathieu: Mtot = 1, G = 1, Etot = -1/4."
  [bs]
  (let [mtot (double (mtot bs))
        bs (map (fn [b] (assoc b :m (/ (:m b) mtot))) bs)]
    (let [eg (double (energy bs))]
      (assert (< eg 0))
      (let [scalefactor (double (/ -0.25 eg))]
        (map (fn [b]
               (let [new-r (v* (/ 1.0 scalefactor) (:r b))
                     new-v (v* (Math/sqrt scalefactor) (:v b))]
                 (assoc b :r new-r :v new-v)))
             bs)))))

(defn- choose-r
  "Returns a radius uniformly distributed in volume within rmax."
  [^Random rng rmax]
  (let [rmax (double rmax)
        rmax2 (* rmax rmax)
        rmax3 (* rmax2 rmax)
        r3 (* rmax3 (.nextDouble rng))]
    (Math/pow r3 (/ 1.0 3.0))))

(defn hot-spherical
  "Returns a system of n equal-mass bodies that has positions and
  velocities uniformly distributed in a sphere."
  [rng n]
  (let [n (int n)
        m (/ (double 1.0) n)
        vmax (Math/sqrt (/ 5.0 6.0))
        rmax (double (/ 6.0 5.0))]
    (letfn [(next-body [id]
              (let [rmag (choose-r rng rmax)
                    vmag (choose-r rng vmax)]
                {:id id :m m :t 0.0 :r (uniform-on-sphere rng rmag) :v (uniform-on-sphere rng vmag)}))]
      (to-standard-units (to-com-frame (map next-body (range n)))))))