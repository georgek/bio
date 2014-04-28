(asdf:defsystem "gk-bio"
  :description "Stuff for bioinformatics."
  :version "0.1"
  :author "George Kettleborough"
  :licence "GNU GPLv3"
  :depends-on ("cl-ppcre" "gzip-stream" "sqlite" "cl-who")
  :serial t
  :components ((:file "packages")
               (:file "utils")
               (:file "generic")
               (:file "bio")
               (:file "genetic-code")
               (:file "sequence")
               (:file "dna")
               (:file "io")
               (:file "db")
               (:file "graphics")
               (:file "strand-bias")))
