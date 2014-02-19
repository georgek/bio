(in-package :gk-bio)

;;; input

(defun char-to-base (char)
  (the (integer 0 3) (logand (ash (char-code char) -1) 3)))

(defun base-to-char (base)
  (aref nucleotides base))

(defun acid-to-char (acid)
  (if (numberp acid)
      (aref amino-acids acid)
      #\*))

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

(defun read-fasta-file (filename)
  "Reads fasta file at FILENAME and outputs list of sequences."
  (let ((sequences (list (make-instance 'dna-sequence))))
    (with-open-file (filein filename)
      (loop for line = (read-line filein nil)
         while line do
           (cond
             ((char= (elt line 0) #\>)
              (when (> (length (bases (car sequences))) 0)
                (push (make-instance 'dna-sequence) sequences))
              (setf (name (car sequences)) (string-trim "> " line)))
             (t
              (loop for char across line do
                   (when (char-dna-basep char)
                    (push-to-sequence (car sequences) (char-to-base char))))))))
    sequences))

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

