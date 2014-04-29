(in-package :gk-bio)

(defun binom (n k)
  (loop for i from 1 to k
     with r = 1
     do
       (setf r (* r (/ (- (1+ n) i) i)))
     finally
       (return r)))

(defun fact (n)
  (loop for i from 1 below n
     for fac = n then (* fac i)
     finally (return fac)))

(defun approx-fact (n)
  (* (sqrt (* 2 pi n))
     (expt (/ n (exp 1)) n)))

(defun approx-log-fact (n &key (more-precision nil))
  "Fast approximation to log(n!) using Stirling's approximation."
  (setf n (coerce n 'float))
  (if (< n 10)
      (log (fact n))
      (let ((lnf (/ (log (* 2.0 pi)) 2.0)))
        (incf lnf (* (log n) (+ n 0.5)))
        (decf lnf n)
        (incf lnf (/ 1.0 (* 12.0 n)))
        (when more-precision
          (decf lnf (/ 1.0 (* 360.0 (expt n 3.0))))
          (incf lnf (/ 1.0 (* 1260.0 (expt n 5.0))))
          (decf lnf (/ 1.0 (* 1280.0 (expt n 7.0)))))
        lnf)))

(defun approx-log-binom (n k)
  (- (approx-log-fact n)
     (+ (approx-log-fact k) (approx-log-fact (- n k)))))

(defun table-probability (a b c d)
  "Prob for table  a b
                   c d"
  (/ (* (binom (+ a b) a)
        (binom (+ c d) c))
     (binom (+ a b c d) (+ a c))))

(defun approx-table-probability (a b c d)
  "Prob for table  a b
                   c d"
  (exp (- (+ (approx-log-binom (+ a b) a)
             (approx-log-binom (+ c d) c))
          (approx-log-binom (+ a b c d) (+ a c)))))

(defun table-probability2 (a b c d)
  "Approximate version of `table-probability'. Not very accurate."
  (let ((r1 (+ a b))
        (r2 (+ c d))
        (c1 (+ a c))
        (c2 (+ b d))
        (n (+ a b c d)))
    (exp (- (+ (- (approx-log-fact r1)
                  (+ (approx-log-fact a) (approx-log-fact b)))
               (- (approx-log-fact r2)
                  (+ (approx-log-fact c) (approx-log-fact d))))
            (- (approx-log-fact n)
               (+ (approx-log-fact c1) (approx-log-fact c2)))))))

(defun fisher-exact-test (a b c d &key (tails :both) (approx nil))
  "Fisher's exact test for table  a b
                                  c d"
  ;; put smallest value in position a
  (flet ((swap-tails ()
           (case tails (:less :greater) (:greater :less) (t :both))))
    (when (or (< b a) (< d c))
      (rotatef a b) (rotatef c d) (setf tails (swap-tails)))
    (when (< c a)
      (rotatef a c) (rotatef b d) (setf tails (swap-tails)))
    (let* ((prob-function (if approx #'table-probability2 #'table-probability))
           (p0 (funcall prob-function a b c d))
           (p1 p0))
      (loop until (zerop a)
         with p = p1
         do
           (incf b) (incf c)
           (setf p (* (/ (* a d) (* b c)) p))
           (incf p1 p)
           (decf a) (decf d))
      (ecase tails
        (:less
         p1)
        (:greater
         (- 1 p1 (- p0)))
        (:both
         (* 2 (min p1 (- 1 p1 (- p0)))))))))

(defun strand-bias (af ar df dr)
  (float (fisher-exact-test af ar df dr :tails :both)))

