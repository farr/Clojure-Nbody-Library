(ns nbody.hermite
  (:use nbody.vec nbody.gravity))

(defrecord Particle [id ^double m ^double t ^double tnext
                     ^doubles r
                     ^doubles v
                     ^doubles a
                     ^doubles j])

(defn map->particle
  "Turns a map with :id, :m, :t, :r, :v members into a particle."
  [m]
  (new Particle
       (:id m)
       (:m m)
       (:t m)
       (:t m)
       (:r m)
       (:v m)
       (double-array 3 0.0)
       (double-array 3 0.0)))

(defn- compare-tnext
  "Compare two particles by tnext."
  [^Particle a ^Particle b]
  (let [ta (double (.tnext a))
        tb (double (.tnext b))]
    (cond
     (< ta tb) -1
     (= ta tb) (compare a b)
     :else 1)))

(defn- predict
  "Fills in rpred and vpred with the predicted position and velocity
  of the given body."
  [^Particle p tpred rpred vpred]
  (let [rpred (doubles rpred)
        vpred (doubles vpred)
        r (doubles (.r p))
        v (doubles (.v p))
        a (doubles (.a p))
        j (doubles (.j p))
        dt (- (double tpred) (double (.t p)))]
    (dotimes [i 3]
      (aset rpred i
            (+ (aget r i)
               (* dt (+ (aget v i)
                        (* (/ dt (double 2.0))
                           (+ (aget a i)
                              (* (/ dt (double 3.0))
                                 (aget j i))))))))
      (aset vpred i
            (+ (aget v i)
               (* dt (+ (aget a i)
                        (* (/ dt (double 2.0))
                           (aget j i)))))))))

(defn- accum-acc-and-jerk
  "Returns [a j] (acceleration and jerk) on b due to the bodies bs."
  [^Particle b bs]
  (let [rpred (double-array 3)
        vpred (double-array 3)
        t (.t b)
        r (.r b)
        v (.v b)
        m (.m b)
        a (double-array 3 0.0)
        j (double-array 3 0.0)
        atot (double-array 3 0.0)
        jtot (double-array 3 0.0)]
    (doseq [^Particle b2 bs]
      (when (not (= b2 b))
        (predict b2 t rpred vpred)
        (acc-and-jerk m r v (.m b2) rpred vpred a j)
        (dotimes [i 3]
          (aset atot i (+ (aget atot i) (aget a i)))
          (aset jtot i (+ (aget jtot i) (aget j i))))))
    [atot jtot]))

(defn- combine-acc-and-jerk
  "Combines [a j] pairs in a sequence into the total [a j]."
  [ajs]
  (reduce
   (fn [[^doubles atot ^doubles jtot] [^doubles a ^doubles j]]
     [(amap atot i aa
            (+ (aget atot i) (aget a i)))
      (amap jtot i jj
            (+ (aget jtot i) (aget j i)))])
   ajs))

(defn- split-particles
  "Partitions a particle sequence into a number of separate sequences
  that can then be processed.  The number of separate sequences is
  always at least twice the number of available processors, to ensure
  good saturation."
  [ps]
  (let [n (count ps)
        nproc (.. Runtime getRuntime availableProcessors)
        n (int (Math/floor (/ n (* 2 nproc))))]
    (partition-all n ps)))

(defn- compute-acc-and-jerk
  "Computes the acceleration and jerk due to bs on b.  Produces the
  correct result even if b is a member of bs.  Tries to saturate the
  available number of processors by computing the acc and jerk in
  parallel."
  [b bs]
  (combine-acc-and-jerk
   (pmap (fn [b2] (accum-acc-and-jerk b b2))
         (split-particles bs))))

(defn- prediction-timescale
  "Returns an estimate of the timescale that would make the predictor
  error equal to epsabs*step-length."
  [^doubles r ^doubles rnew ^doubles rpred dt epsrel]
  (let [step-distance (distance r rnew)
        pred-distance (distance rnew rpred)]
    (if (= pred-distance 0.0)
      (* dt (Math/pow 2.0 (/ 1.0 3.0)))
      (* dt (Math/pow (/ (* epsrel step-distance)
                         pred-distance)
                      (/ 1.0 3.0))))))

(defn- advance-particle
  "Returns a new particle that represents the time-advancement of b
  due to the effects of the bs.  epsrel is used to set the timestep,
  such that the prediction error in the step is epsrel*(change in r
  over step)."
  [^Particle b bs epsrel]
  (let [tnext (double (.tnext b))
        t (double (.t b))
        dt (double (- tnext t))
        r (doubles (.r b))
        v (doubles (.v b))
        a (doubles (.a b))
        j (doubles (.j b))
        rp (double-array 3 0.0)
        vp (double-array 3 0.0)]
    (predict b tnext rp vp)
    (let [[^doubles anew ^doubles jnew]
          (compute-acc-and-jerk (assoc b :r rp :t tnext :v vp) bs)]
      (let [vnew (doubles
                  (amap v i vnew
                        (+ (aget v i)
                           (* dt
                              (+ (* (double 0.5)
                                    (+ (aget a i) (aget anew i)))
                                 (* dt
                                    (* (/ (double 1.0) (double 12.0))
                                       (- (aget j i) (aget jnew i)))))))))
            rnew (amap r i rnew
                       (+ (aget r i)
                          (* dt
                             (+ (* (double 0.5)
                                   (+ (aget v i) (aget vnew i)))
                                (* dt
                                   (* (/ (double 1.0) (double 12.0))
                                      (- (aget a i) (aget anew i))))))))
            dt (prediction-timescale r rnew rp dt epsrel)]
        (assoc b :t tnext :tnext (+ tnext dt) :r rnew :v vnew :a anew :j jnew)))))

(defn- first-step-timescale
  [^Particle b epsabs]
  (let [j (.j b)
        rpred (double-array 3 0.0)
        vpred (double-array 3 0.0)]
    (predict b (+ (.t b) 1.0) rpred vpred)
    (let [step-distance (distance (.r b) rpred)
          pred-distance (/ (norm j) 6.0)]
      (Math/sqrt (/ (* epsabs step-distance) pred-distance)))))

(defn- setup-before-integration
  [bs epsrel]
  (map
   (fn [^Particle b]
     (let [[a j]
           (compute-acc-and-jerk b bs)]
       (let [dt (first-step-timescale (assoc b :a a :j j) epsrel)]
         (assoc b :a a :j j :tnext (+ (.t b) dt)))))
   bs))

(defn advance
  "Advances a set of bodies from their current time to tstop (which
  must be in the future).  epsrel is used to set the individual body
  timesteps."
  [bs tstop epsrel]
  (let [bs (setup-before-integration (map map->particle bs) epsrel)]
    (let [new-bs
          (loop [bs (apply sorted-set-by compare-tnext bs)]
            (let [^Particle b (first bs)
                  other-bs (disj bs b)]
              (if (> (.tnext b) tstop)
                bs
                (let [new-b (advance-particle b other-bs epsrel)]
                  (recur (conj other-bs new-b))))))]
      (reduce (fn [done-bs b]
                (cons (advance-particle (assoc b :tnext tstop) (disj new-bs b) epsrel)
                      done-bs))
              ()
              new-bs))))