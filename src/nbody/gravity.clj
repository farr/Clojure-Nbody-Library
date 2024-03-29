(ns nbody.gravity
  (:use nbody.vec))

(definline acc-and-jerk!
  "Fills a and j with the acceleration and jerk on body 1 due to body 2."
  [r1 v1 m2 r2 v2 a j]
  `(let [r1# (doubles ~r1)
         v1# (doubles ~v1)
         m2# (double ~m2)
         r2# (doubles ~r2)
         v2# (doubles ~v2)
         a# (doubles ~a)
         j# (doubles ~j)]
     (let [rdv# (double
                 (areduce r1# i# rdv# (double 0.0)
                          (+ rdv#
                             (* (- (aget r2# i#)
                                   (aget r1# i#))
                                (- (aget v2# i#)
                                   (aget v1# i#))))))]
       (let [r# (distance r1# r2#)
             r22# (* r# r#)
             r3# (* r# r22#)
             r5# (* r3# r22#)]
         (dotimes [i# 3]
           (aset a# i#
                 (/ (* m2# (- (aget r2# i#) (aget r1# i#)))
                    r3#))
           (aset j# i#
                  (* m2#
                     (- (/ (- (aget v2# i#) (aget v1# i#))
                           r3#)
                        (* (double 3.0)
                           (* rdv#
                              (/ (- (aget r2# i#) (aget r1# i#))
                                 r5#)))))))))))

(definline acc!
  "Fills a with the acceleration on body 1 due to body 2."
  [r1 m2 r2 a]
  `(let [r1# (doubles ~r1)
         m2# (double ~m2)
         r2# (doubles ~r2)
         a# (doubles ~a)]
     (let [r# (double (distance r1# r2#))
           r22# (* r# r#)
           r3# (* r22# r#)]
       (dotimes [i# 3]
         (aset a# i# (/ (* m2# (- (aget r2# i#) (aget r1# i#)))
                        r3#))))))