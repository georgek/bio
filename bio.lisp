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
  (aref amino-acids acid))

