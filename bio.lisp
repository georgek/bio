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

(defun char-to-base (char)
  (the (integer 0 3) (logand (ash (char-code char) -1) 3)))

(defun base-to-char (base)
  (aref nucleotides base))

(defun base-complement (base)
  (the (integer 0 3) (logxor 2 base)))

(defun read-seq (stream &optional name)
  (let ((seq (make-instance 'dna-sequence :name name)))
    (loop for char = (read-char stream nil)
       while char do
         (when (find char "acgt")
           (vector-push-extend (char-to-base char) (bases seq))))
    seq))

(defun read-seq-file (filename &optional name)
  (with-open-file (filein filename)
    (read-seq filein name)))

(defparameter bases-per-chunk 10)
(defparameter chunks-per-line 6)

(defun write-seq (seq &key (pretty nil) (stream t))
  (loop for base across (bases seq)
     for i from 0 do
       (when (and pretty (= 0 (mod i (* bases-per-chunk chunks-per-line))))
         (format stream "~9d" (1+ i)))
       (when (and pretty (= 0 (mod i bases-per-chunk)))
         (format stream " "))
       (write-char (base-to-char base) stream)
       (when (and pretty (= 0 (mod (1+ i) (* bases-per-chunk chunks-per-line))))
         (format stream "~%"))))

(defun write-seq-file (seq filename &key (pretty nil))
  (with-open-file (fileout filename :direction :output)
    (write-seq seq :pretty pretty :stream fileout)
    (fresh-line fileout)))

(defmacro other (form o1 o2)
  `(if (eq ,form ,o1) ,o2 ,o1))

(defun seq-complement (seq)
  (make-instance 'dna-sequence
                 :name (name seq)
                 :direction (other (direction seq) :5p :3p)
                 :bases (map 'vector #'base-complement (bases seq))))

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

(defun triplet-to-acid (base1 base2 base3)
  (aref gene-code base1 base2 base3))

(defun acid-to-char (acid)
  (if (numberp acid)
      (aref amino-acids acid)
      #\*))

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
    (sort (reduce #'nconc orfs) #'< :key #'car)))

