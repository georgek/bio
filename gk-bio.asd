(asdf:defsystem "gk-bio"
  :description "Stuff for bioinformatics."
  :version "0.1"
  :author "George Kettleborough"
  :licence "GNU GPLv3"
  :serial t
  :components ((:file "packages")
               (:file "dna")
               (:file "io")
               (:file "bio")))
