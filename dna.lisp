(in-package :gk-bio)

(defclass dna-sequence ()
  ((name
    :initarg :name
    :initform nil
    :accessor name
    :documentation "Name of sequence.")
   (bases
    :initarg :bases
    :initform (make-array 100
                          :element-type '(integer 0 3)
                          :adjustable t
                          :fill-pointer 0)
    :accessor bases)
   (direction
    :initarg :direction
    :initform :5p
    :accessor direction
    :documentation "The direction of the sequence (:5p or :3p).")))

(defmethod print-object ((object dna-sequence) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "with ~A bases" (length (bases object)))))

(defmethod length ((sequence dna-sequence))
  (length (bases sequence)))

(defmethod elt ((sequence dna-sequence) index)
  (aref (bases sequence) index))

(defmethod push-to-sequence ((sequence dna-sequence) base)
  (assert (typep base '(integer 0 3)))
  (vector-push-extend base (bases sequence)))

(defun char-dna-basep (char)
  (find char "acgtACGT"))

(defun base-complement (base)
  (the (integer 0 3) (logxor 2 base)))

(defmacro other (form o1 o2)
  `(if (eq ,form ,o1) ,o2 ,o1))

(defun seq-complement (seq)
  (make-instance 'dna-sequence
                 :name (name seq)
                 :direction (other (direction seq) :5p :3p)
                 :bases (map 'vector #'base-complement (bases seq))))

(defun triplet-to-acid (base1 base2 base3)
  (aref gene-code base1 base2 base3))

(defun find-orfs (dna-seq min-length)
  (let ((bases (bases dna-seq))
        (lengths (make-array 6 :initial-element 0))
        (orfs (map 'vector (lambda (n) (list (cons n 0))) #(0 1 2 0 1 2))))
    (loop for i from 0
       for f = (mod i 3)
       for b = (+ f 3)
       while (< i (- (length (bases dna-seq)) 2)) do
         ;; forward
         (if (triplet-to-acid (aref bases i) (aref bases (+ i 1))
                              (aref bases (+ i 2)))
             ;; not a stop codon
             (incf (aref lengths f) 3)
             (progn                     ; stop codon
               (when (>= (aref lengths f) min-length)
                 (setf (cdr (car (aref orfs f))) (aref lengths f))
                 (push (cons i 0) (aref orfs f)))
               (setf (car (aref orfs f)) (cons (1+ i) 0))
               (setf (aref lengths f) 0)))
         ;; backward
         (if (triplet-to-acid (base-complement (aref bases (+ i 2)))
                              (base-complement (aref bases (+ i 1)))
                              (base-complement (aref bases i)))
             ;; not a stop codon
             (incf (aref lengths b) 3)
             (progn                     ; stop codon
               (when (>= (aref lengths b) min-length)
                 (setf (cdr (car (aref orfs b))) (- (aref lengths b)))
                 (push (cons i 0) (aref orfs b)))
               (setf (car (aref orfs b)) (cons (1+ i) 0))
               (setf (aref lengths b) 0))))
    (setf orfs (map 'vector (lambda (f) (if (< (abs (cdar f)) min-length)
                                            (cdr f)
                                            f))
                    orfs))
    (sort (reduce #'nconc orfs) #'> :key (lambda (orf) (abs (cdr orf))))))

