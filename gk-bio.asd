(asdf:defsystem "gk-bio"
  :description "Stuff for bioinformatics."
  :version "0.1"
  :author "George Kettleborough"
  :licence "GNU GPLv3"
  :serial t
  :components ((:file "packages")
               (:file "utils")
               (:file "generic")
               (:file "bio")
               (:file "genetic-code")
               (:file "dna")
               (:file "io")))
