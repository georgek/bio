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
                            db "select A,C,G,T
                              from pileupnd where animal = ? and day = ?
                              order by chromosome, position;"
                            animal-id day)
            do
              (push (cons (format nil "~A-d~D" animal-name day)
                          sequence)
                    sequences)))
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

