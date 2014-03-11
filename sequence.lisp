(in-package :gk-bio)

;;; generic sequence which can contain any character
(defclass seq ()
  ((name
    :initarg :name
    :initform nil
    :accessor name
    :documentation "Name of sequence.")
   (characters
    :initarg :elements
    :initform (make-array 100
                          :element-type 'character
                          :adjustable t
                          :fill-pointer 0)
    :accessor characters
    :documentation "The elements of the sequence.")
   (direction
    :initarg :direction
    :initform :5p
    :accessor direction
    :documentation "The direction of the sequence (:5p or :3p).")))

(defmethod initialize-instance :after ((seq seq) &key)
  (loop for i from 0 below (length seq) do
       (setf (aref (characters seq) i)
             (char-downcase (aref (characters seq) i)))))

(defmethod print-object ((object seq) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "with ~A characters" (length (characters object)))))

(defmethod length ((sequence seq))
  (length (characters sequence)))

(defmethod elt ((sequence seq) index)
  (aref (characters sequence) index))

(defmethod push-to-sequence ((sequence seq) character)
  (vector-push-extend (char-downcase character) (characters sequence)))

;;; protein, a sequence of amino acids (and nothing else)
(defclass protein-sequence ()
  ((name
    :initarg :name
    :initform nil
    :accessor name
    :documentation "Name of protein.")
   (acids
    :initarg :acids
    :initform (make-array 100
                          :element-type '(integer 0 20)
                          :adjustable t
                          :fill-pointer 0)
    :accessor acids)))

(defmethod length ((sequence protein-sequence))
  (length (acids sequence)))

(defmethod elt ((sequence protein-sequence) index)
  (aref (acids sequence) index))

(defmethod push-to-sequence ((sequence protein-sequence) base)
  (assert (typep base '(integer 0 20)))
  (vector-push-extend base (acids sequence)))
