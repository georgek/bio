(in-package :gk-bio)


(defclass protein-sequence ()
  ((name
    :initarg :name
    :initform nil
    :accessor name
    :documentation "Name of protein.")
   (acids
    :initarg :acids
    :initform (make-array 100
                          :element-type '(integer 0 19)
                          :adjustable t
                          :fill-pointer 0)
    :accessor acids)))

