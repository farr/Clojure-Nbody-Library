(ns nbody.gravity
  (:use nbody.vec))

(definline acc-and-jerk
  "Fills a and j with the acceleration and jerk on body 1 due to body 2."
  [m1 r1 v1 m2 r2 v2 a j]
  `(let [m1# (double ~m1)
         r1# (doubles ~r1)
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
         (dotimes [i# (int 3)]
           (aset a# i#
                 (/ (* m1# (* m2# (- (aget r2# i#) (aget r1# i#))))
                    r3#))
           (aset j# i#
                  (* (* m1# m2#)
                      (- (/ (- (aget v2# i#) (aget v1# i#))
                             r3#)
                         (* (double 3.0)
                             (* rdv#
                                (/ (- (aget r2# i#) (aget r1# i#))
                                   r5#)))))))))))