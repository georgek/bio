(in-package :gk-bio)

(defgeneric length (sequence))

(defmethod length (sequence)
  (cl:length sequence))

(defgeneric elt (sequence index))

(defmethod elt (sequence index)
  (cl:elt sequence index))

(defgeneric push-to-sequence (sequence base)
  (:documentation "Adds BASE to the end of SEQUENCE."))

