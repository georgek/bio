(in-package :gk-bio)

(defconstant nucleotides #(#\a #\c #\t #\g))
(defconstant an 0)
(defconstant cn 1)
(defconstant gn 3)
(defconstant tn 2)

(defconstant amino-acids #(#\a #\r #\n #\d #\c #\e #\q #\g #\h #\i
                           #\l #\k #\m #\f #\p #\s #\t #\w #\y #\v #\*))

(defun gene-table ()
  (let ((gt (make-array '(4 4 4))))
   (loop for acid across
        "ffllssssyy**cc*wllllpppphhqqrrrriiimttttnnkkssrrvvvvaaaaddeegggg"
      for one across
        "ttttttttttttttttccccccccccccccccaaaaaaaaaaaaaaaagggggggggggggggg"
      for two across
        "ttttccccaaaaggggttttccccaaaaggggttttccccaaaaggggttttccccaaaagggg"
      for three across
        "tcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcag"
      do
        (setf (aref gt
                    (position one nucleotides)
                    (position two nucleotides)
                    (position three nucleotides))
              (position acid amino-acids)))
   gt))

(defconstant gene-code
  #3A(((11 2 2 11) (16 16 16 16) (9 9 9 12) (1 15 15 1))
      ((6 8 8 6) (14 14 14 14) (10 10 10 10) (1 1 1 1))
      ((20 18 18 20) (15 15 15 15) (10 13 13 10) (20 4 4 17))
      ((5 3 3 5) (0 0 0 0) (19 19 19 19) (7 7 7 7))))

