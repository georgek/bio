(in-package :gk-bio)

;;; input

(defun char-to-base (char)
  (the (integer 0 3) (logand (ash (char-code char) -1) 3)))

(defun base-to-char (base &key (upper nil))
  (let ((char (aref nucleotides base)))
    (if upper
        (char-upcase char)
        char)))

(defun acid-to-char (acid)
  (aref amino-acids acid))

(defun read-seq (stream &optional name)
  (let ((seq (make-instance 'dna-sequence :name name)))
    (loop for char = (read-char stream nil)
       while char do
         (when (char-dna-basep char)
           (push-to-sequence seq (char-to-base char))))
    seq))

(defun seq (string &optional name)
  (read-seq (make-string-input-stream string) name))

(defun read-seq-file (filename &optional name)
  (with-open-file (filein filename)
    (read-seq filein name)))

(defun file-is-gzip (filename)
  "Determines if a file is a GZIP file based on its magic number."
  (with-open-file (filein filename :element-type 'unsigned-byte)
    (and (= (read-byte filein) #x1f)
         (= (read-byte filein) #x8b))))

(defmacro with-open-file-maybe-gzip ((stream filespec &key (direction :input))
                                     &body body)
  `(if (file-is-gzip ,filespec)
       ,(macroexpand `(with-open-gzip-file (,stream ,filespec :direction ,direction)
                        ,@body))
       ,(macroexpand `(with-open-file (,stream ,filespec :direction ,direction)
                        ,@body))))

(defun read-fasta-file (filename)
  "Reads fasta file at FILENAME and outputs list of sequences."
  (let ((sequences (list (make-instance 'seq))))
    (with-open-file-maybe-gzip (filein filename)
      (loop for line = (read-line filein nil)
         while line do
           (when (> (length line) 0)
             (cond
               ((char= (elt line 0) #\>)
                (when (> (length (characters (car sequences))) 0)
                  (push (make-instance 'seq) sequences))
                (setf (name (car sequences)) (string-trim "> " line)))
               (t
                (loop for char across line do
                     (push-to-sequence (car sequences) char)))))))
    (if (> (length (characters (car sequences))) 0)
        (nreverse sequences)
        (nreverse (cdr sequences)))))

(defun file-read-error (filename line-number)
  (error (format nil "Error in file ~S at line ~D~%" filename line-number)))

(defun read-fastq-file (filename)
  "Reads fastq file at FILENAME and outputs a list of sequences."
  (let ((sequences (list)))
    (with-open-file-maybe-gzip (filein filename)
      (loop with line-number = 1
         for line = (read-line filein nil)
         while line do
           (push (make-instance 'dna-sequence) sequences)
           (if (char= (elt line 0) #\@)
               (setf (name (car sequences)) (string-trim "@ " line))
               (file-read-error filename line-number))
           (setf line (read-line filein nil))
           (incf line-number)
           (loop for char across line do
                (if (char-dna-basep char)
                    (push-to-sequence (car sequences) (char-to-base char))
                    (file-read-error filename line-number)))
           (setf line (read-line filein nil))
           (incf line-number)
           (unless (char= (elt line 0) #\+)
             (file-read-error filename line-number))
           (setf line (read-line filein nil))
           (incf line-number)
           ;; this line has the quality scores
           (incf line-number)))
    (nreverse sequences)))

;;; output

(defparameter items-per-chunk 10)
(defparameter chunks-per-line 6)

(defun pretty-print-seq (seq &key (stream t) (key #'identity))
  "Prints a vector in a pretty way to STREAM using function KEY to transform
each element."
  (loop for item across seq
     for i from 0 do
       (when (= 0 (mod i (* items-per-chunk chunks-per-line)))
         (format stream "~9d" (1+ i)))
       (when (= 0 (mod i items-per-chunk))
         (format stream " "))
       (write-char (funcall key item) stream)
       (when (= 0 (mod (1+ i) (* items-per-chunk chunks-per-line)))
         (format stream "~%"))))

(defgeneric write-seq (seq &key pretty stream)
  (:documentation "Write SEQ object to STREAM.  If PRETTY is true it is done
  in a pretty way.")  )

(defmethod write-seq ((seq dna-sequence) &key (pretty nil) (stream t))
  (if pretty
      (pretty-print-seq (bases seq) :stream stream :key #'base-to-char)
      (loop for item across (bases seq) do
           (write-char (base-to-char item) stream))))

(defmethod write-seq ((seq protein-sequence) &key (pretty nil) (stream t))
  (if pretty
      (pretty-print-seq (acids seq) :stream stream :key #'acid-to-char)
      (loop for item across (acids seq) do
           (write-char (acid-to-char item) stream))))

(defun write-seq-file (seq filename &key (pretty nil))
  (with-open-file (fileout filename :direction :output)
    (write-seq seq :pretty pretty :stream fileout)
    (fresh-line fileout)))

(defun write-fasta-file (sequences stream)
  (when (atom sequences)
    (setf sequences (list sequences)))
  (loop for sequence in sequences do
       (format stream ">~A" (name sequence))
       (loop for base across (bases sequence)
          for i from 0 do
            (when (= (mod i 80) 0)
              (write-char #\Newline stream))
            (write-char (base-to-char base :upper t) stream))
       (fresh-line stream)))
