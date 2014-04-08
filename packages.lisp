(in-package :cl-user)

(defpackage :gk-bio
  (:use :common-lisp :gzip-stream)
  (:shadow :length :elt))
