(in-package :gk-bio)

(defparameter name-length 12)

(defun html-char (character)
  (case character
    (#\- "&#8209;")
    (t (string character))))

(defun gen-html (sequences stream &optional (reference nil))
  (if (null reference) (setf reference (consensus sequences)))
  (who:with-html-output (stream nil :prologue t)
    (:html
     (:head (:title "Test")
            (:link :rel "stylesheet" :href "test.css"))
     (:body
      (let ((len (min (reduce #'min (mapcar #'length sequences))
                      (length reference))))
        (loop for sequence in sequences do
             (who:htm
              (:div
               (:span :class "name"
                      (who:fmt (format nil "~~~DA" name-length)
                               (subseq (name sequence) 0
                                       (min name-length
                                            (length sequence)))))
               (loop for i from 0 below len
                  for character = (elt sequence i) do
                    (who:htm
                     (:span :class
                            (if (eq character (elt reference i))
                                character
                                "snp")
                            (who:fmt "~A" (html-char character))))))))))))
  (fresh-line stream)
  (values))

