(in-package :gk-bio)

(defgeneric length (sequence))

(defmethod length (sequence)
  (cl:length sequence))

(defgeneric elt (sequence index))

(defmethod elt (sequence index)
  (cl:elt sequence index))

