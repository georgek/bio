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
         (when (find char "acgt")
           (vector-push-extend (char-to-base char) (bases seq))))
    seq))

(defun seq (string &optional name)
  (read-seq (make-string-input-stream string) name))

(defun read-seq-file (filename &optional name)
  (with-open-file (filein filename)
    (read-seq filein name)))

;;; output

(defparameter items-per-chunk 10)
(defparameter chunks-per-line 6)

(defun pretty-print-seq (seq &key (stream t) (key #'identity))
  (loop for item across seq
     for i from 0 do
       (when (= 0 (mod i (* items-per-chunk chunks-per-line)))
         (format stream "~9d" (1+ i)))
       (when (= 0 (mod i items-per-chunk))
         (format stream " "))
       (write-char (funcall key item) stream)
       (when (= 0 (mod (1+ i) (* items-per-chunk chunks-per-line)))
         (format stream "~%"))))

(defgeneric write-seq (seq &key pretty stream))

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
