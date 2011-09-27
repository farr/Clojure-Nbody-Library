(ns nbody.hermite
  (:use nbody.vec nbody.gravity))

(definline compare-double-array
  [a1 a2]
  `(let [a1# (doubles ~a1)
         a2# (doubles ~a2)]
     (let [n1# (int (alength a1#))
           n2# (int (alength a2#))]
       (let [c# (int (compare n1# n2#))]
         (if (= c# 0)
           (loop [i# (int 0)]
             (cond
              (>= i# n1#) 0
              (= (aget a1# i#) (aget a2# i#)) (recur (+ i# 1))
              (< (aget a1# i#) (aget a2# i#)) -1
              :else 1))
           c#)))))

(defrecord Particle [id ^double m ^double t ^double tnext
                     ^doubles r ^doubles v
                     ^doubles a ^doubles j]
  java.lang.Comparable
  (compareTo [this other]
    (let [^Particle this this
          ^Particle other other]
      (let [c (int (compare (.tnext this) (.tnext other)))]
        (if (= c (int 0))
          (let [c (int (compare (.m this) (.m other)))]
            (if (= c (int 0))
              (let [c (int (compare (.t this) (.t other)))]
                (if (= c (int 0))
                  (let [c (int (compare-double-array (.r this) (.r other)))]
                    (if (= c (int 0))
                      (let [c (int (compare-double-array (.v this) (.v other)))]
                        (if (= c (int 0))
                          (let [c (int (compare-double-array (.a this) (.a other)))]
                            (if (= c (int 0))
                              (let [c (int (compare-double-array (.j this) (.j other)))]
                                c)
                              c))
                          c))
                      c))
                  c))
              c))
          c)))))

(defn map->particle
  "Takes a map with components :id, :m, :t, :r, :v, and converts it to
  a Particle."
  [m]
  (if (instance? Particle m)
    m
    (new Particle
         (:id m) (:m m) (:t m) (:t m)
         (:r m) (:v m) (double-array 3 0.0) (double-array 3 0.0))))

(defn print-particle
  "Outputs a particle in PSDF format."
  [^Particle p]
  (io!
   (print "--- !Particle\n")
   (print "id: " (.id p) "\n")
   (print "m: " (.m p) "\n")
   (print "t: " (.t p) "\n")
   (print "t_max: " (.tnext p) "\n")
   (let [r (doubles (.r p))] (print "r: [" (aget r 0) ", " (aget r 1) ", " (aget r 2) "]\n"))
   (let [v (doubles (.v p))] (print "v: [" (aget v 0) ", " (aget v 1) ", " (aget v 2) "]\n"))
   (let [a (doubles (.a p))] (print "a: [" (aget a 0) ", " (aget a 1) ", " (aget a 2) "]\n"))
   (let [j (doubles (.j p))] (print "j: [" (aget j 0) ", " (aget j 1) ", " (aget j 2) "]\n"))))

(defn predict!
  "Fills rp and vp with the predicted position of p at time tp."
  [^Particle p tp rp vp]
  (let [r (doubles (.r p))
        v (doubles (.v p))
        a (doubles (.a p))
        j (doubles (.j p))
        dt (- (double tp) (double (.t p)))
        tnext (double (.tnext p))
        rp (doubles rp)
        vp (doubles vp)]
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

(definline vec=
  "Returns true only if the given double-arrays are equal."
  [v w]
  `(let [v# (doubles ~v)
         w# (doubles ~w)
         n# (int (alength v#))]
     (and (= n# (int (alength w#)))
          (loop [i# (int 0)]
            (or (>= i# n#)
                (let [vi# (aget v# i#)
                      wi# (aget w# i#)]
                  (and (= vi# wi#)
                       (recur (+ i# (int 1))))))))))

(defn total-acc-and-jerk
  "Returns the total acc and jerk on p due to ps, computed at (:t b),
  predicting the positions of ps as necessary."
  [^Particle p ps]
  (let [atot (doubles (double-array 3 0.0))
        jtot (doubles (double-array 3 0.0))
        a (doubles (double-array 3 0.0))
        j (doubles (double-array 3 0.0))
        r (doubles (.r p))
        v (doubles (.v p))
        t (double (.t p))
        rp (doubles (double-array 3 0.0))
        vp (doubles (double-array 3 0.0))]
    (doseq [^Particle p2 ps]
      (let [m2 (double (.m p2))]
        (predict! p2 t rp vp)
        (when (not (vec= r rp))
          (acc-and-jerk! r v m2 rp vp a j)
          (dotimes [i 3]
            (aset atot i (+ (aget atot i) (aget a i)))
            (aset jtot i (+ (aget jtot i) (aget j i)))))))
    [atot jtot]))

(defn acc-jerk-timescale
  "Returns a timescale based on the given acc and jerk."
  [a j sf]
  (* (* 2.0 (Math/sqrt (double sf))) (/ (norm a) (norm j))))

(defn prediction-failure-timescale
  "Returns the timescale that would be estimated from the failure of
  the predicted velocity to match the corrected velocity."
  [vp vnew a h sf]
  (let [dv (double (distance vp vnew))
        delta-v (double (* h (norm a)))
        sf (double sf)]
    (if (= dv 0.0)
      (* h 2.01) ; Double stepsize if predictor is too accurate. 
      (* h (Math/sqrt (* sf (/ delta-v dv)))))))

(defn advanced-particle
  "Returns a particle that is b advanced to (.tnext b) in the system
  of bs.  sf is a timestep safety factor: the next timestep size will
  be chosen by sf*timescale."
  [^Particle b bs sf]
  (let [rp (doubles (double-array 3 0.0))
        vp (doubles (double-array 3 0.0))
        r (doubles (.r b))
        v (doubles (.v b))
        a (doubles (.a b))
        j (doubles (.j b))
        tp (double (.tnext b))
        t (double (.t b))
        h (- tp t)]
    (predict! b tp rp vp)
    (let [bp (assoc b :t tp :r rp :v vp)
          [^doubles ap ^doubles jp] (total-acc-and-jerk bp bs)]
      (let [vnew (doubles (amap v i vnew
                                (+ (aget v i)
                                   (* h
                                      (+ (* (double 0.5) (+ (aget a i) (aget ap i)))
                                         (* h
                                            (* (/ (double 1.0) (double 12.0))
                                               (- (aget j i) (aget jp i)))))))))]
        (let [rnew (amap r i rnew
                         (+ (aget r i)
                            (* h
                               (+ (* (double 0.5) (+ (aget v i) (aget vnew i)))
                                  (* h
                                     (* (/ (double 1.0) (double 12.0))
                                        (- (aget a i) (aget ap i))))))))
              hnew (prediction-failure-timescale vp vnew a h sf)]
          (assoc bp :t tp :r rnew :v vnew :a ap :j jp :tnext (+ tp hnew)))))))

(defn setup
  "Sets up a particle array for integration.  sf is the timestep
  safety factor."
  [ps sf]
  (let [ps (map map->particle ps)]
    (map
     (fn [p]
       (if (= (:t p) (:tnext p))
         (let [[a j] (total-acc-and-jerk p ps)
               t (:t p)
               h (acc-jerk-timescale a j sf)
               tnext (+ t h)]
           (assoc p :a a :j j :tnext tnext))
         p))
     ps)))

(defn particle-can-advance?
  "True iff :tnext is smaller than the given time."
  [^Particle p t]
  (let [tn (double (.tnext p))
        t (double t)]
    (> tn t)))

(defn minimum-tnext
  [ps]
  (reduce (fn [mt ^Particle p]
            (let [mt (double mt)
                  pt (double (.tnext p))]
              (min mt pt)))
          Double/POSITIVE_INFINITY
          ps))

(defn shared-timestep-advance
  "Advance all bodies by the same timestep to the given time."
  [ps tstop sf]
  (let [ps (setup ps sf)]
    (loop [ps ps]
      (let [tn (minimum-tnext ps)]
        (if (>= tn tstop)
          (pmap #(advanced-particle (assoc % :tnext tstop) ps sf) ps)
          (recur (pmap #(advanced-particle (assoc % :tnext tn) ps sf) ps)))))))

(defn internal-advance
  "Advances ps beginning at tstart to tstop."
  [ps tstart tstop tup sf]
  (if-let [these-ps (seq (take-while (fn [p] (and (particle-can-advance? p tstop)
                                                  (not (particle-can-advance? p tup))))
                                     ps))]
    (let [new-ps (pmap #(advanced-particle (assoc % :tnext tstop) ps sf) these-ps)]
      (reduce (fn [ps [old-p new-p]]
                (conj (disj ps old-p) new-p))
              ps (map (fn [x y] [x y]) these-ps new-ps)))
    (let [thalf (+ tstart (/ (- tstop tstart) 2.0))
          half-ps (internal-advance ps tstart thalf tstop sf)]
      (recur half-ps thalf tstop tup sf))))

(defn advance
  "Advances bs using individual timesteps to tstop, using sf as the
  timestep safety factor."
  [ps tstop sf]
  (let [ps (setup ps sf)
        tstart (:t (first ps))
        ps (apply sorted-set ps)]
    (internal-advance ps tstart tstop Double/POSITIVE_INFINITY sf)))