(in-package :gk-bio)

;;; for square matrices
(defclass matrix ()
  ((names
    :initarg :names
    :initform nil
    :accessor names
    :documentation "List of names of columns.")
   (values
    :initarg :vals
    :initform nil
    :accessor vals
    :documentation "The matrix.")))

(defmethod print-object ((object matrix) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "size ~A" (length (names object)))))

(defun make-matrix (names)
  (make-instance 'matrix :names names
                 :vals (make-array `(,(length names) ,(length names))
                                   :initial-element nil)))

(defun copy-matrix (matrix)
  (make-instance 'matrix :vals (vals matrix) :names (names matrix)))

(defun matrix-elt (matrix i j &key (test #'eql))
  "Get element by name of column and row."
  (aref (vals matrix)
        (position i (names matrix) :test test)
        (position j (names matrix) :test test)))

(defun (setf matrix-elt) (val matrix i j &key (test #'eql))
  (setf (aref (vals matrix)
              (position i (names matrix) :test test)
              (position j (names matrix) :test test))
        val))

(defun print-matrix (matrix stream &key (delimiter #\,))
  (loop for cname on (names matrix) do
       (format stream "~A" (car cname))
       (unless (endp (cdr cname))
         (format stream "~C" delimiter)))
  (fresh-line stream)
  (let ((ilim (1- (array-dimension (vals matrix) 0)))
        (jlim (1- (array-dimension (vals matrix) 1))))
   (loop for i from 0 to ilim do
        (loop for j from 0 to jlim
           for val = (aref (vals matrix) i j) do
             (if val
                 (format stream "~A" val)
                 (format stream "-"))
             (unless (= j jlim)
               (format stream "~C" delimiter)))
        (fresh-line stream))))

(defun manhattan (a b)
  "A and B lists of numbers.  Calculates Manhattan distance between A and B."
  (reduce #'+ (mapcar (lambda (a b) (abs (- a b))) a b)))

(defun position-variation-distance (pos1 pos2)
  (let* ((sum1 (float (reduce #'+ pos1)))
         (sum2 (float (reduce #'+ pos2)))
         (norm1 (mapcar (lambda (b) (if (zerop sum1) 0.0 (/ b sum1))) pos1))
         (norm2 (mapcar (lambda (b) (if (zerop sum2) 0.0 (/ b sum2))) pos2)))
    (/ (manhattan norm1 norm2)
       2)))

(defun sequence-variation-distance (seq1 seq2)
  "SEQ1 and SEQ2 must be lists of equal length of lists of numbers of equal
  length. The distance is the sum of half the Manhattan distances between
  corresponding positions in the sequences.  Is equivalent to the Manhattan
  distances if only one base per position, ie. each sublist has one nonzero."
  (reduce #'+ (mapcar #'position-variation-distance seq1 seq2)))

(defun sequence-variation-matrix (db)
  "Makes sequence variation matrix for every sample in db."
  (let ((sequences (list))
        distance-matrix)
    (loop for (animal-id animal-name)
       in (sqlite:execute-to-list
           db "select id,name from animals;")
       do
         (loop for (day) in (sqlite:execute-to-list
                             db "select distinct day from pileup
                                 where animal = ?;"
                             animal-id)
            for sequence = (sqlite:execute-to-list
                            db "select Af,Ar,Cf,Cr,Gf,Gr,Tf,Tr
                                from pileup where animal = ? and day = ?
                                order by chromosome, position;"
                            animal-id day)
            do
              (dbg :seqvarmat "Removing bias...~%")
              (setf sequence (mapcar (lambda (p) (apply #'filter-bias p))
                                     sequence))
              (setf sequence (mapcar (lambda (p) (apply #'double-to-single-strand p))
                                     sequence))
              (push (cons (format nil "~A-d~D" animal-name day)
                          sequence)
                    sequences)))
    (setf sequences (nreverse sequences))
    (setf distance-matrix (make-matrix (mapcar #'car sequences)))
    (loop for name in (names distance-matrix) do
         (setf (matrix-elt distance-matrix name name :test #'equal) 0))
    (loop for seqc on sequences
       for seq1 = (first seqc)
       do
         (dbg :seqvarmat "~A~%" (car seq1))
         (loop for seq2 in (rest seqc)
            for distance = (sequence-variation-distance
                            (cdr seq1)
                            (cdr seq2))
            do
              (setf (matrix-elt distance-matrix
                                (car seq1) (car seq2)
                                :test #'equal)
                    distance)
              (setf (matrix-elt distance-matrix
                                (car seq2) (car seq1)
                                :test #'equal)
                    distance)))
    distance-matrix))

(defun consensus-base (counts)
  "Returns base char with highest count.  Counts should be a list of a,c,g,t
counts."
  (let ((alist (mapcar #'cons counts (list #\a #\c #\g #\t))))
    (cdr (assoc (reduce #'max counts) alist :test #'=))))

(defun window-var (seq1 seq2 window-size filter-function)
  (flet ((pospr (pos) (apply #'double-to-single-strand
                             (apply filter-function pos))))
    (loop with length = (min (length seq1) (length seq2))
       for i from 0 below length by window-size
       collect
         (sequence-variation-distance
          (mapcar #'pospr (subseq seq1 i (min (+ i window-size) length)))
          (mapcar #'pospr (subseq seq2 i (min (+ i window-size) length)))))))

(defun print-win-var-diff (stream seq1 seq2 window-size filter1 filter2)
  "For comparing filtering methods."
  (let ((windows1 (window-var seq1 seq2 window-size filter1))
        (windows2 (window-var seq1 seq2 window-size filter2)))
    (loop with length = (min (length seq1) (length seq2))
       with diffs = (mapcar #'- windows1 windows2)
       for diff in diffs
       for i from 0 below length by window-size
       do
         (format stream "~5D - ~5D: ~7,3F~%"
                 (1+ i) (min (+ i window-size) length) diff))))

(defun print-win-var-seq-diff (stream refseq seq1 seq2 window-size filter)
  "For comparing sequences to a reference."
  (let ((windows1 (window-var refseq seq1 window-size filter))
        (windows2 (window-var refseq seq2 window-size filter)))
    (loop with length = (min (length refseq) (length seq1) (length seq2))
       with diffs = (mapcar #'- windows1 windows2)
       for diff in diffs
       for i from 0 below length by window-size
       do
         (format stream "~5D - ~5D: ~7,3F~%"
                 (1+ i) (min (+ i window-size) length) diff))))
