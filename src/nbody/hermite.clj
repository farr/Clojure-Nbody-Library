(ns nbody.hermite
  (:use nbody.vec nbody.gravity))

(defrecord Particle [id ^double m ^double t ^double tnext
                     ^doubles r ^doubles v
                     ^doubles a ^doubles j])

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
  [a j]
  (/ (norm a) (norm j)))

(defn acc-jerk-predicted-timescale
  "Returns a timescale based on the "
  [a j ap jp h]
  (let [a (doubles a)
        j (doubles j)
        ap (doubles ap)
        jp (doubles jp)
        h (double h)
        ; The following is v_corr - v_pred, and is O(h^3), despite appearances.
        dv (double (norm (amap a i dv
                               (- (* (double 0.5)
                                     (* h (- (aget ap i) (aget a i))))
                                  (* (/ (double 1.0) (double 12.0))
                                     (* (* h h)
                                        (+ (* (double 5.0) (aget j i))
                                           (aget jp i))))))))]
    (if (= dv 0.0)
      (* h 2.1) ; Increase h enough to move up in block-power-of-two sizes.
      (Math/sqrt (/ (* (double (norm a)) h)
                    dv)))))

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
              hnew (* sf (acc-jerk-predicted-timescale a j ap jp h))]
          (assoc bp :t tp :r rnew :v vnew :a ap :j jp :tnext (+ tp hnew)))))))

(defn setup
  "Sets up a particle array for integration.  sf is the timestep
  safety factor.  Ensures that no particle's timestep takes it
  beyond tstop."
  [ps sf]
  (let [ps (map map->particle ps)]
    (map
     (fn [p]
       (let [[a j] (total-acc-and-jerk p ps)
             t (:t p)
             h (* sf (acc-jerk-timescale a j))
             tnext (+ t h)]
         (assoc p :a a :j j :tnext tnext)))
     ps)))

(defn compare-by-tnext
  "Compares particles first by their :tnext value, and then by the standard comparison ordering."
  [^Particle p1 ^Particle p2]
  (let [t1 (double (.tnext p1))
        t2 (double (.tnext p2))]
    (let [c (int (compare t1 t2))]
      (if (= c 0)
        (compare p1 p2)
        c))))