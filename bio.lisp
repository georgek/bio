(in-package :gk-bio)


(defclass protein-sequence ()
  ((name
    :initarg :name
    :initform nil
    :accessor name
    :documentation "Name of protein.")
   (acids
    :initarg :acids
    :initform (make-array 100
                          :element-type '(integer 0 19)
                          :adjustable t
                          :fill-pointer 0)
    :accessor acids)))

(defun median (sequence &key (order #'<))
  (let ((sequence (sort sequence order)))
    (if (evenp (length sequence))
        (/ (+ (elt sequence (/ (length sequence) 2))
              (elt sequence (1- (/ (length sequence) 2)))) 2)
        (elt sequence (floor (/ (length sequence) 2))))))

(defun n50 (sequences)
  "Calculates n50 statistic for a list of sequences."
  (let* ((lengths (sort (mapcar #'length sequences) #'>))
         (total (reduce #'+ lengths)))
    (dbg :n50 "~A~%" lengths)
    (loop with sum = 0
       for length in lengths
       do (incf sum length)
       while (< sum (/ total 2))
       finally (return length))))

(defun gc-content (sequences)
  (let ((gc-counts (mapcar (lambda (n-counts) (+ (nth (char-to-base #\c) n-counts)
                                                 (nth (char-to-base #\g) n-counts)))
                           (mapcar #'nucleotide-count sequences))))
    (/ (reduce #'+ gc-counts)
       (reduce #'+ (mapcar #'length sequences)))))

