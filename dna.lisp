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
    (format stream "~A with ~A bases" (name object) (length (bases object)))))

(defmethod seq-to-dna (seq)
  (let ((dna (make-instance 'dna-sequence
                            :name (name seq)
                            :bases (make-array (length seq)
                                               :element-type '(integer 0 3)
                                               :adjustable t
                                               :fill-pointer (length seq)))))
    (loop for i from 0 below (length seq) do
         (if (char-dna-basep (aref (characters seq) i))
             (setf (aref (bases dna) i) (char-to-base (aref (characters seq) i)))
             (error "~C is not a nucleotide.~%" (aref (characters seq) i))))
    dna))

(defmethod length ((sequence dna-sequence))
  (length (bases sequence)))

(defmethod elt ((sequence dna-sequence) index)
  (aref (bases sequence) index))

(defmethod push-to-sequence ((sequence dna-sequence) base)
  (check-type base (integer 0 3))
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
                 :bases (nreverse (map 'vector #'base-complement (bases seq)))))

(defun triplet-to-acid (base1 base2 base3)
  (aref gene-code base1 base2 base3))

(defun seq-translate (seq &optional (frame 0))
  (let ((protein (make-instance 'protein-sequence)))
    (loop for i from frame to (- (length seq) 3) by 3
       for acid = (triplet-to-acid (elt seq i)
                                   (elt seq (+ i 1))
                                   (elt seq (+ i 2)))
       do
         (dbg :trans "~A ~A ~A => ~A~%" (elt seq i) (elt seq (+ i 1)) (elt seq (+ i 2))
              acid)
         (push-to-sequence protein acid))
    protein))

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

(defun nucleotide-count (dna-sequence)
  "Returns a list of the four nucleotide counts.  Use char-to-base to index by
character."
  (let ((counts (make-list 4 :initial-element 0)))
    (loop for base across (bases dna-sequence) do
         (incf (nth base counts)))
    counts))

