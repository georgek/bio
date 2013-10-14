(in-package :gk-bio)

(defun lmax (list)
  (reduce #'max list))

(defun plmax (list)
  (position (lmax list) list))

(defun align (seq1 seq2
              &key (match-score 5) (mismatch-score -2) (space-score -6))
  (let* ((m (length seq1))
         (n (length seq2))
         (matrix (make-array (list m n) :initial-element nil))
         (seq1a (make-array (max m n) :element-type 'character
                            :initial-element #\-))
         (seq2a (make-array (max m n) :element-type 'character
                            :initial-element #\-)))
    (labels ((sigma (a b) (if (= a b) match-score mismatch-score))
             (scores (i j)
               (cond
                 ((= i 0) (list (* j space-score)))
                 ((= j 0) (list (* i space-score)))
                 (t (aref matrix (1- i) (1- j))))))
      ;; fill in matrix
      (loop for i below m do
           (loop for j below n do
                (setf (aref matrix i j)
                      (list (+ (lmax (scores i j))
                               (sigma (elt seq1 i) (elt seq2 j)))
                            (+ (lmax (scores i (1+ j)))
                               space-score)
                            (+ (lmax (scores (1+ i) j))
                               space-score)))))
      ;; find alignment
      (loop with i = (1- m)
         with j = (1- n)
         for k from (1- (max m n)) downto 0
         while (>= i 0) while (>= j 0) do
           (format t "i: ~a, j: ~a, k: ~a~%" i j k)
           (case (plmax (aref matrix i j))
             (0 (setf (aref seq1a k) (base-to-char (elt seq1 i)))
                (setf (aref seq2a k) (base-to-char (elt seq2 j)))
                (decf i)
                (decf j))
             (1 (setf (aref seq1a k) (base-to-char (elt seq1 i)))
                (setf (aref seq2a k) #\-)
                (decf i))
             (2 (setf (aref seq1a k) #\-)
                (setf (aref seq2a k) (base-to-char (elt seq2 j)))
                (decf j))))
      (format t "Sequence 1: ~s~%Sequence 2: ~s~%" seq1a seq2a))))

