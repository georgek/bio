(in-package :gk-bio)

(defgeneric length (sequence))

(defmethod length (sequence)
  (cl:length sequence))

(defgeneric elt (sequence index))

(defmethod elt (sequence index)
  (cl:elt sequence index))

(defgeneric push-to-sequence (sequence object)
  (:documentation "Adds OBJECT to the end of SEQUENCE."))

(defgeneric get-char (sequence index)
  (:documentation "Returns the printable character corresponding to this
  position."))

