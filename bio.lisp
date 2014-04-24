(in-package :gk-bio)



(defun median (sequence &key (order #'<))
  (let ((sequence (sort sequence order)))
    (if (evenp (length sequence))
        (/ (+ (elt sequence (/ (length sequence) 2))
              (elt sequence (1- (/ (length sequence) 2)))) 2)
        (elt sequence (floor (/ (length sequence) 2))))))

(defun mode (list)
  (let ((counts (make-hash-table)))
    (loop for item in list do
         (if (gethash item counts)
             (incf (gethash item counts))
             (setf (gethash item counts) 1)))
    (loop with max-length = 0 with max
       for item being the hash-keys in counts do
         (when (> (gethash item counts) max-length)
           (setf max item
                 max-length (gethash item counts)))
       finally (return max))))

(defun n50 (sequences)
  "Calculates n50 statistic for a list of sequences."
  (let* ((lengths (sort (mapcar #'length sequences) #'>))
         (total (reduce #'+ lengths)))
    (dbg :n50 "~A~%" lengths)
    (loop with sum = 0
       for length in lengths
       do (incf sum length)
       while (< sum (/ total 2))
       finally (return length))))

(defun gc-content (sequences)
  (let ((gc-counts (mapcar (lambda (n-counts) (+ (nth (char-to-base #\c) n-counts)
                                                 (nth (char-to-base #\g) n-counts)))
                           (mapcar #'nucleotide-count sequences))))
    (/ (reduce #'+ gc-counts)
       (reduce #'+ (mapcar #'length sequences)))))

(defun consensus (sequences)
  "Builds a consensus sequence from the list of sequences.  The consensus will
  be the length of the shortest sequence and each position is determined by a
  majority vote."
  (let* ((length (reduce #'min (mapcar #'length sequences)))
         (consensus (make-instance 'seq)))
    (loop for i from 0 below length do
         (push-to-sequence consensus (mode (mapcar (lambda (s) (elt s i)) sequences))))
    consensus))

(defun all-equal (collection &key (test #'eq))
  (every #'identity
         (map 'list (lambda (e) (funcall test e (elt collection 0)))
              collection)))

(defun find-repeats (sequences n &optional (accepted-bases "acgt"))
  (setf accepted-bases (string-downcase accepted-bases))
  (loop for seq in sequences
     for bases = (characters seq)
     for length = (length seq)
     do
       (loop for i from 0 upto (- length (1+ (* 2 n))) do
            (when (and
                   (position (char-downcase (elt bases i)) accepted-bases)
                   (not (eq (elt bases i) (elt bases (+ i n))))
                   (not (eq (elt bases (+ i n)) (elt bases (+ i n 1))))
                   (eq (elt bases i) (elt bases (+ i n 1)))
                   (all-equal (subseq bases i (+ i n)))
                   (all-equal (subseq bases (+ i n 1) (+ i n 1 n))))
              (format t "Seq: ~A, pos: ~D, found: ~{~A~}~A~{~A~}~%"
                      (name seq) (+ i n)
                      (map 'list #'identity (subseq bases i (+ i n)))
                      (elt bases (+ i n))
                      (map 'list #'identity (subseq bases (+ i n 1) (+ i n 1 n))))))))
