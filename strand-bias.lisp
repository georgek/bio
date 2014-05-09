(in-package :gk-bio)

(defun binom (n k)
  (dbg :binom "binom ~A ~A~%" n k)
  (loop for i from 1 to k
     with r = 1
     do
       (setf r (* r (/ (- (1+ n) i) i)))
     finally
       (return r)))

(defun fact (n)
  (loop for fac = (max 1 n) then (* fac i)
     for i from 1 below n
     finally (return fac)))

(defun approx-fact (n)
  (* (sqrt (* 2 pi n))
     (expt (/ n (exp 1)) n)))

(defun approx-log-fact (n &key (more-precision nil))
  "Fast approximation to log(n!) using Stirling's approximation."
  (setf n (coerce n 'long-float))
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
  (dbg :tab-prob "Table a: ~A, b: ~A, c: ~A, d: ~A~%" a b c d)
  (/ (* (binom (+ a b) a)
        (binom (+ c d) c))
     (binom (+ a b c d) (+ a c))))

(defun approx-table-probability (a b c d)
  "Prob for table  a b
                   c d"
  (exp (- (+ (approx-log-binom (+ a b) a)
             (approx-log-binom (+ c d) c))
          (approx-log-binom (+ a b c d) (+ a c)))))

(defmacro with-coercion (things type &body body)
  `(let
     ,(mapcar (lambda (thing) `(,thing (coerce ,thing ,type)))
               things)
     ,@body))

(defun table-probability2 (a b c d)
  "Approximate version of `table-probability'. Not very accurate."
  (with-coercion (a b c d) 'long-float
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
                (+ (approx-log-fact c1) (approx-log-fact c2))))))))

(defun rotate-table (a b c d)
  ;; rotate to make A the smallest
  (loop while (/= a (min a b c d)) do
       (rotatef a b d c))
  (values a b c d))

(defun transpose-table (a b c d)
  ;; transpose to make args to binom smaller
  (when (< b c)
    (rotatef b c))
  (values a b c d))

(defun other-extreme (a b c d)
  "Returns the most extreme table in the other tail"
  (let ((swap-rows nil)
        (swap-cols nil))
    (when (< (+ c d) (+ a b))
      (rotatef a c) (rotatef b d)
      (setf swap-rows t))
    (when (< (+ b d) (+ a c))
      (rotatef a b) (rotatef c d)
      (setf swap-cols t))
    ;; now the small row and col totals are at the top and left
    (let ((diff (- a (min b c))))
      (decf a diff) (decf d diff)
      (incf b diff) (incf c diff))
    (when swap-rows
      (rotatef a c) (rotatef b d))
    (when swap-cols
      (rotatef a b) (rotatef c d))
    (multiple-value-call #'transpose-table
      (rotate-table a b c d))))

;;; method from Jerrold H. Zar
;;; "A fast and efficient algorithm for the Fisher exact test".
;;; doi: 10.3758/BF03202590
(defun fisher-exact-test (a b c d &key (tails :both) (approx nil))
  "Fisher's exact test for table  a b
                                  c d"
  ;; put smallest value in position a
  (flet ((swap-tails ()
           (case tails (:less :greater) (:greater :less) (t :both))))
    (multiple-value-bind (a b c d)
        (rotate-table a b c d)
      ;; rotating and transposing doesn't affect the table-probabilty, but
      ;; transposing changes the tail
      (when (< b c)
        (rotatef b c)
        (setf tails (swap-tails)))
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
          (:greater
           p1)
          (:less
           (- 1 (- p1 p0)))
          (:both
           (multiple-value-bind (a b c d)
               (other-extreme a b c d)
             (let* ((p (funcall prob-function a b c d))
                    (p2 0))
               (loop while (and (< p p0) (>= b 0) (>= c 0)) do
                    (incf p2 p)
                    (incf a) (incf d)
                    (setf p (* (/ (* b c) (* a d)) p))
                    (decf b) (decf c))
               (+ p1 p2)))))))))

(defun strand-bias (af ar df dr)
  (- 1 (float (fisher-exact-test af ar df dr :tails :both))))

