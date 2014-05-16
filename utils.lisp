(in-package :gk-bio)

(defvar *dbg-ids* nil
  "Identifiers used by DBG")

(defun dbg (id format-string &rest args)
  "Print debugging info if (DEBUG ID) has been specified"
  (when (member id *dbg-ids*)
    (fresh-line *debug-io*)
    (apply #'format *debug-io* format-string args)))

(defun dbgo (object)
  (format *debug-io* "~A~%" object)
  object)

(defun set-debug (&rest ids)
  "Start dbg output on the given ids"
  (setf *dbg-ids* (union ids *dbg-ids*)))

(defun undebug (&rest ids)
  "Stop dbg on ids. With no ids, stop dbg altogether"
  (setf *dbg-ids* (if (null ids) nil
                      (set-difference *dbg-ids* ids))))

(defun dbg-on-p (id)
  (member id *dbg-ids*))

(defun sanitise-filename (filename)
  (string-trim "_-." (ppcre:regex-replace-all
                      "[^a-zA-Z0-9_\\-.]+" filename "_")))

(defun head (list &optional (n 10))
  (loop repeat n
     for element in list collect element))

(defun compose (&rest funs)
  (setf funs (nreverse funs))
  (lambda (&rest args)
    (loop with val = (apply (first funs) args)
       for fun in (rest funs) do
         (setf val (funcall fun val))
       finally
         (return val))))

(defun range (beg end &optional (step 1))
  "Returns range of numbers between beg and end."
  (assert (and (<= beg end) (> step 0)))
  (loop for i from beg to end by step collecting i))

(defun split-string (string delimiter &key (omit-nulls t))
  (assert (stringp string))
  (assert (characterp delimiter))
  (let ((splits (list)))
   (loop for pos = (position delimiter string)
      while pos do
        (push (subseq string 0 pos) splits)
        (setf string (subseq string (1+ pos)))
      finally (push string splits))
   (when omit-nulls
     (setf splits (delete "" splits :test #'equal)))
   (nreverse splits)))

