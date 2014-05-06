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
     (:head (:title "Haplotypes")
            (:link :rel "stylesheet" :href "test.css"))
     (:body
      (:script :src "wz_tooltip/wz_tooltip.js" :type "text/javascript")
      (let ((len (min (reduce #'min (mapcar #'length sequences))
                      (length reference))))
        (loop for sequence in sequences do
             (who:htm
              (:div
               (:span :class "name"
                      :onmouseover (format nil "Tip(\"~A\");" (name sequence))
                      :onmouseout "UnTip();"
                      (who:fmt (format nil "~~~DA" name-length)
                               (subseq (name sequence) 0
                                       (min name-length
                                            (length (name sequence))))))
               (loop for i from 0 below len
                  for character = (elt sequence i) do
                    (who:htm
                     (:span :class
                            (if (eq character (elt reference i))
                                character
                                "snp")
                            :onmouseover (format nil "Tip(~A);" (1+ i))
                            :onmouseout "UnTip();"
                            (who:fmt "~A" (html-char character))))))))))))
  (fresh-line stream)
  (values))

(defparameter nuc-colours #("green" "blue" "red" "yellow"))

(defun tikz-mult-algn (dna-sequences stream)
  "Makes TikZ code to show sequence logo of aligned sequences."
  (format stream "\\begin{tikzpicture}[font=\\sffamily\\bfseries\\Huge,x=1pt,y=1pt]~%")
  (let ((n (length dna-sequences))
        (length (reduce #'min dna-sequences :key #'length)))
    (loop for i from 0 below length
       for x from 0 by 15
       for counts = (nucleotide-count
                     (make-instance 'dna-sequence :bases
                                    (map 'vector
                                         (lambda (seq) (aref (bases seq) i))
                                         dna-sequences)))
       do
         (loop for j from 0 to 3
            for prop = (/ (float (aref counts j)) n)
            with y = 0
            do
              (unless (zerop prop)
               (format stream "\\node[yscale=~,3F,anchor=base,color=~A] at (~D,~,3Fex) {~C};~%"
                       prop
                       (aref nuc-colours j)
                       x
                       y
                       (base-to-char j :upper t))
               (incf y (* 4.0 prop))))))
    (format stream "\\end{tikzpicture}")
    (fresh-line stream))
