(in-package :gk-bio)

(defun binom (n k)
  (loop for i from 1 to k
     with r = 1
     do
       (setf r (* r (/ (- (1+ n) i) i)))
     finally
       (return r)))

(defun table-probability (a b c d)
  "Prob for table  a b
                   c d"
  (/ (* (binom (+ a b) a)
        (binom (+ c d) c))
     (binom (+ a b c d) (+ a c))))

(defun fisher-exact-test (a b c d &key (tails :both))
  "Fisher's exact test for table  a b
                                  c d"
  ;; put smallest value in position a
  (flet ((swap-tails (tail)
           (case tails (:less :greater) (:greater :less) (t :both))))
    (when (or (< b a) (< d c))
      (rotatef a b) (rotatef c d) (setf tails (swap-tails tails)))
    (when (< c a)
      (rotatef a c) (rotatef b d) (setf tails (swap-tails tails)))
    (let* ((p0 (table-probability a b c d))
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

