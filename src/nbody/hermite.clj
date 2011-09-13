(ns nbody.hermite
  (:use nbody.vec nbody.gravity))

(defrecord Particle [id ^double m ^double t ^double tnext
                     ^doubles r ^doubles v
                     ^doubles a ^doubles j])

(defn map->particle
  "Convert a map with id, m, t, r, v entries to a particle."
  [b]
  (if (instance? Particle b)
    b
    (new Particle
         (:id b) (:m b) (:t b) (:t b) (:r b) (:v b)
         (double-array 3 0.0) (double-array 3 0.0))))

(defn predict!
  "Fills the given rp vp arrays with the predicted position of p at
  time tp."
  [^Particle p tp rp vp]
  (let [t (double (.t p))
        tp (double tp)
        dt (- tp t)
        rp (doubles rp)
        vp (doubles vp)
        r (doubles (.r p))
        v (doubles (.v p))
        a (doubles (.a p))
        j (doubles (.j p))]
    (dotimes [i 3]
      (aset rp i
            (+ (aget r i)
               (* dt
                  (+ (aget v i)
                     (* (/ dt (double 2.0))
                        (+ (aget a i)
                           (* (/ dt (double 3.0))
                              (aget j i))))))))
      (aset vp i
            (+ (aget v i)
               (* dt
                  (+ (aget a i)
                     (* (/ dt (double 2.0))
                        (aget j i)))))))))

(definline arr=
  [x y]
  `(let [x# (doubles ~x)
         y# (doubles ~y)]
     (loop [i# (int 0)]
       (if (>= i# (alength x#))
         true
         (let [xx# (aget x# i#)
               yy# (aget y# i#)]
           (and (= xx# yy#)
                (recur (+ i# (int 1)))))))))

(defn total-acc-and-jerk
  "Returns the acc and jerk on b due to the bs.  bs can contain b."
  [^Particle b bs]
  (let [r (doubles (.r b))
        t (double (.t b))
        v (doubles (.v b))
        acc (double-array 3 0.0)
        jerk (double-array 3 0.0)
        rp (double-array 3 0.0)
        vp (double-array 3 0.0)
        a (double-array 3 0.0)
        j (double-array 3 0.0)]
    (doseq [^Particle b2 bs]
      (let [m2 (double (.m b2))]
        (predict! b2 t rp vp)
        (when (not (arr= rp r))
          (acc-and-jerk! r v m2 rp vp a j)
          (dotimes [i 3]
          (aset acc i (+ (aget acc i) (aget a i)))
          (aset jerk i (+ (aget jerk i) (aget j i)))))))
    [acc jerk]))

(defn predictor-error-timescale
  "Returns the timescale that would have the predictor error for a
  particle be a fraction epsrel of the total distance covered in a
  timestep."
  [rold r rp dt epsrel]
  (let [dr (distance r rold)
        d (distance r rp)]
    (if (= d 0.0)
      (* dt (Math/pow 2.0 0.25))
      (let [new-dt (* dt (Math/pow (double (/ (* epsrel dr) d)) (/ 1.0 3.0)))]
        new-dt))))

(defn first-step-timescale
  "Returns a predictor error timescale for the first step of a particle."
  [^doubles v ^doubles a ^doubles j epsrel]
  (let [dr (norm
            (amap v i ret (+ (aget v i) (+ (* (double 0.5) (aget a i)) (/ (aget j i) (double 6.0))))))
        d (/ (norm j) 6.0)]
    (let [dt (Math/sqrt (double (/ (* epsrel dr) d)))]
      dt)))

(defn advance-particle
  "Returns an advanced particle."
  [^Particle b bs epsrel]
  (let [r (doubles (.r b))
        v (doubles (.v b))
        a (doubles (.a b))
        j (doubles (.j b))
        rp (doubles (double-array 3 0.0))
        vp (doubles (double-array 3 0.0))
        t (double (.t b))
        tp (double (.tnext b))
        dt (double (- tp t))]
    (predict! b tp rp vp)
    (let [[^doubles ap ^doubles jp]
          (total-acc-and-jerk (assoc b :t tp :r rp :v vp) bs)]
      (let [vnew (doubles (amap v i ret
                                (+ (aget v i)
                                    (* dt
                                        (+ (* (double 0.5) (+ (aget a i) (aget ap i)))
                                           (* dt
                                              (* (/ (double 1.0) (double 12.0))
                                                 (- (aget j i) (aget jp i)))))))))
            rnew (doubles (amap r i ret
                                (+ (aget r i)
                                   (* dt
                                      (+ (* (double 0.5) (+ (aget vnew i) (aget v i)))
                                         (* dt
                                            (* (/ (double 1.0) (double 12.0))
                                               (- (aget a i) (aget ap i)))))))))]
        (assoc b :t tp :r rnew :v vnew :a ap :j jp :tnext (+ tp (predictor-error-timescale r rnew rp dt epsrel)))))))

(defn compare-tnext
  "Comparison function that compares two Particles based on their
  tnext values."
  [^Particle p1 ^Particle p2]
  (let [t1 (double (.tnext p1))
        t2 (double (.tnext p2))]
    (let [c (int (compare t1 t2))]
      (if (= c 0)
        (compare p1 p2)
        c))))

(defn setup-for-advance
  "Preps bodies for advance by the integrator."
  [bs epsrel]
  (let [bs (map map->particle bs)]
    (map
     (fn [b]
       (let [[a j] (total-acc-and-jerk b bs)]
         (assoc b :a a :j j :tnext (+ (:t b) (first-step-timescale (:v b) a j epsrel)))))
     bs)))

(defn advance
  "Returns an advanced system of bodies, stopping at exactly tstop."
  [bs tstop epsrel]
  (let [bs (apply sorted-set-by compare-tnext (setup-for-advance bs epsrel))]
    (let [new-bs
          (loop [bs bs]
            (let [^Particle b (first bs)]
              (if (> (.tnext b) tstop)
                bs
                (let [rem (disj bs b)]
                  (recur (conj rem (advance-particle b rem epsrel)))))))]
      (map #(advance-particle (assoc % :tnext tstop) new-bs epsrel) new-bs))))