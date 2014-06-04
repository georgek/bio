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

(defun hamming (seq1 seq2)
  "SEQ1 and SEQ2 should be sequences."
  (reduce #'+ (map 'list (lambda (a b) (if (equal a b) 0 1)) seq1 seq2)))

(defun consensus-base (counts)
  "Returns base char with highest count.  Counts should be a list of a,c,g,t
counts."
  (let ((alist (mapcar #'cons counts (list #\a #\c #\g #\t))))
    (cdr (assoc (reduce #'max counts) alist :test #'=))))

(defun discretise-position (counts)
  "Returns the position but with the maximum count set to 1 and others set to
0."
  (let ((max (reduce #'max counts)))
    (mapcar (lambda (c) (if (= c max) 1 0)) counts)))

(defun consensus-variation-matrix (db)
  (let ((consensuses (list))
        distance-matrix)
    (loop for (animal-id animal-name)
       in (sqlite:execute-to-list
           db "select id,name from animals;")
       do
         (loop for (day) in (sqlite:execute-to-list
                             db "select distinct day from pileup
                                 where animal = ?;"
                             animal-id)
            do
              (push (cons (format nil "~A-d~D" animal-name day)
                          (get-consensus-sequences2 db animal-id day))
                    consensuses)))
    (setf consensuses (nreverse consensuses))
    (setf distance-matrix (make-matrix (mapcar #'car consensuses)))
    (loop for name in (names distance-matrix) do
         (setf (matrix-elt distance-matrix name name :test #'equal) 0))
    (loop for consensusc on consensuses
       for consensus1 = (first consensusc)
       do
         (loop for consensus2 in (rest consensusc)
              for distance = (reduce #'+ (mapcar
                                          #'hamming
                                          (mapcar #'characters (cdr consensus1))
                                          (mapcar #'characters (cdr consensus2))))
            do
              (setf (matrix-elt distance-matrix
                                (car consensus1) (car consensus2)
                                :test #'equal)
                    distance)
              (setf (matrix-elt distance-matrix
                                (car consensus2) (car consensus1)
                                :test #'equal)
                    distance)))
    distance-matrix))

(defun manhattan (a b)
  "A and B lists of numbers.  Calculates Manhattan distance between A and B."
  (reduce #'+ (mapcar (lambda (a b) (abs (- a b))) a b)))

(defun euclidean-sq (a b)
  (reduce #'+ (mapcar (lambda (a b) (expt (- a b) 2)) a b)))

(defun euclidean (a b)
  (sqrt (euclidean-sq a b)))

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

(defmacro defsimpleclass (name direct-superclasses slots)
  (flet ((add-convenience (slot)
           (nconc slot
                  (list :initarg (intern (symbol-name (car slot)) :keyword)
                        :accessor (intern (concatenate
                                           'string
                                           (symbol-name name) "-"
                                           (symbol-name (car slot))))
                        :initform (case (nth 2 slot)
                                    (integer 0)
                                    (single-float 0.0))))))
    `(defclass ,name ,direct-superclasses
       ,(mapcar #'add-convenience slots))))

(defsimpleclass vphaser-snp ()
  ((chromosome :type string)
   (position :type integer)
   (variant :type (integer 0 3))
   (reference :type (integer 0 3))
   (sb-pval :type single-float)
   (occurrence :type single-float)
   (Af :type integer)
   (Ar :type integer)
   (Cf :type integer)
   (Cr :type integer)
   (Gf :type integer)
   (Gr :type integer)
   (Tf :type integer)
   (Tr :type integer)))

(defmethod print-object ((object vphaser-snp) stream)
  (with-slots (chromosome position Af Ar Cf Cr Gf Gr Tf Tr)
      object
    (print-unreadable-object (object stream :type t)
      (format stream "~@[Chr: ~A ~]Pos: ~D A:~D,~D C:~D,~D G:~D,~D T:~D,~D"
              chromosome position Af Ar Cf Cr Gf Gr Tf Tr))))

(defun read-vphaser-file (filename &optional (chromosome nil))
  "Returns list of vphaser variants. Ignores the LPs at the moment."
  (let ((variants (list)))
    (with-open-file (filein filename)
      (loop for line = (read-line filein nil)
         while line do
           (unless (char= (aref line 0) #\#) ;comment line
             (let ((split (split-string line #\Tab)))
               (cond
                 ((string= (nth 4 split) "snp")
                  (let ((counts (make-array '(8) :initial-element 0)))
                    (loop for count in (nthcdr 6 split)
                       for countsplit = (split-string count #\:)
                       do
                         (case (aref (nth 0 countsplit) 0)
                           (#\A
                            (setf (aref counts 0) (read-from-string
                                                   (nth 1 countsplit)))
                            (setf (aref counts 1) (read-from-string
                                                   (nth 2 countsplit))))
                           (#\C
                            (setf (aref counts 2) (read-from-string
                                                   (nth 1 countsplit)))
                            (setf (aref counts 3) (read-from-string
                                                   (nth 2 countsplit))))
                           (#\G
                            (setf (aref counts 4) (read-from-string
                                                   (nth 1 countsplit)))
                            (setf (aref counts 5) (read-from-string
                                                   (nth 2 countsplit))))
                           (#\T
                            (setf (aref counts 6) (read-from-string
                                                   (nth 1 countsplit)))
                            (setf (aref counts 7) (read-from-string
                                                   (nth 2 countsplit))))))
                    (push (make-instance
                           'vphaser-snp
                           :chromosome chromosome
                           :position (read-from-string (nth 0 split))
                           :variant (char-to-base (aref (nth 1 split) 0))
                           :reference (char-to-base (aref (nth 2 split) 0))
                           :sb-pval (read-from-string (nth 3 split))
                           :occurrence (read-from-string (nth 5 split))
                           :Af (aref counts 0) :Ar (aref counts 1)
                           :Cf (aref counts 2) :Cr (aref counts 3)
                           :Gf (aref counts 4) :Gr (aref counts 5)
                           :Tf (aref counts 6) :Tr (aref counts 7))
                          variants))))))))
    (nreverse variants)))

(defun vphaser-variation-distance (snps1 snps2)
  "Measures distance between sets of SNPs."
  (let ((snp-sites (make-hash-table :test #'equal)))
    ;; We assume that there is only one variant per site so use just the
    ;; variant occurrence rate from vphaser 2. If only one site contains the
    ;; variant then the position variation distance is equal to that
    ;; occurrence rate. Otherwise it is equal to the absolute difference of
    ;; both occurrence rates.
    (loop for snp in (append snps1 snps2)
       for chr = (vphaser-snp-chromosome snp)
       for pos = (vphaser-snp-position snp)
       for occ = (/ (vphaser-snp-occurrence snp) 100)
       do
         (if (gethash (cons chr pos) snp-sites)
             (setf (gethash (cons chr pos) snp-sites)
                   (abs (- (gethash (cons chr pos) snp-sites)
                           occ)))
             (setf (gethash (cons chr pos) snp-sites)
                   occ)))
    (loop for occurrence being the hash-values in snp-sites
       sum occurrence)))

(defun vphaser-variaton-matrix (directories)
  "Calculates the variation matrix using only the variants outputted by
vphaser 2. DIRECTORIES should be a list of strings containing pathnames which
the vphaser outputs can be located, one directory per sample."
  (let* ((sample-names (mapcar (lambda (dir) (car (last (split-string dir #\/))))
                               directories))
         (sample-snps
          (loop for dir in directories collect
               (loop for filename in (directory (format nil "~A/*.fdr.var.txt" dir))
                  for chromosome = (nth 0 (split-string (pathname-name filename) #\.))
                  nconc
                    (read-vphaser-file filename chromosome))))
         (matrix (make-matrix sample-names)))
    (dbg :vphaser-variation-matrix "~A~%" sample-names)
    (loop for name in sample-names do
         (setf (matrix-elt matrix name name :test #'equal) 0))
    (loop for i from 0
       for ic on sample-snps
       do
         (loop for j from (1+ i)
            for jc on (cdr ic)
            for dist = (vphaser-variation-distance (car ic) (car jc))
            do
              (setf (aref (vals matrix) i j) dist)
              (setf (aref (vals matrix) j i) dist)))
    matrix))

