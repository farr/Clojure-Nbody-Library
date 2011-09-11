(ns nbody.vec
  (:import java.util.Random))

(definline dot
  [x y]
  `(let [x# (doubles ~x)
         y# (doubles ~y)]
     (areduce x# i# ret# (double 0.0) (+ ret# (* (aget x# i#)
                                                 (aget y# i#))))))

(definline norm
  [x]
  `(let [x# (doubles ~x)]
     (Math/sqrt (dot x# x#))))

(definline distance-squared
  [x y]
  `(let [x# (doubles ~x)
         y# (doubles ~y)]
     (areduce x# i# ret# (double 0.0)
              (let [dx# (- (aget x# i#) (aget y# i#))]
                (+ ret# (* dx# dx#))))))

(definline distance
  [x y]
  `(Math/sqrt (distance-squared ~x ~y)))

(definline v-
  [x y]
  `(let [x# (doubles ~x)
         y# (doubles ~y)]
     (amap x# i# ret#
           (- (aget x# i#)
              (aget y# i#)))))

(definline v+
  [x y]
  `(let [x# (doubles ~x)
         y# (doubles ~y)]
     (amap x# i# ret#
           (+ (aget x# i#)
              (aget y# i#)))))

(definline v*
  [s v]
  `(let [s# (double ~s)
         v# (doubles ~v)]
     (amap v# i# ret#
           (* s# (aget v# i#)))))

(defn uniform-on-sphere
  "Returns a three-vector uniformly distributed on the sphere of radius r."
  [^Random rng r]
  (let [phi (* 2.0 Math/PI (.nextDouble rng))
        cos-theta (+ -1.0 (* 2.0 (.nextDouble rng)))
        sin-theta (Math/sqrt (- 1.0 (* cos-theta cos-theta)))]
    (double-array
     [(* r (Math/cos phi) sin-theta)
      (* r (Math/sin phi) sin-theta)
      (* r cos-theta)])))