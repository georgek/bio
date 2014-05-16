(in-package :gk-bio)

(defun get-consensus-sequences (db animal-id day)
  (sqlite:execute-non-query
   db
   "create temp table consensus as
  select position, case max(sum(A), sum(C), sum(T), sum(G))
          WHEN sum(A) then 'A'
          WHEN sum(C) then 'C'
          WHEN sum(G) then 'G'
          WHEN sum(T) then 'T'
          END nuc, chromosome
  from pileupnd
  where animal = ? and day = ?
  group by chromosome, position;"
   animal-id day)

  (let (chromosomes
        consensus
        (consensuses (list)))
    (setf chromosomes
          (sqlite:execute-to-list db "select id,name from chromosomes;"))
    (loop for (chr-id chr-name) in chromosomes do
         (setf consensus (mapcar #'car
                                 (sqlite:execute-to-list
                                  db
                                  "select nuc from consensus where chromosome = ?
                                   order by position asc;"
                                  chr-id)))
         (push (make-instance 'seq :name chr-name
                              :characters (map 'vector #'character consensus))
               consensuses))
    (sqlite:execute-non-query db "drop table consensus;")
    (nreverse consensuses)))

(defun get-consensus-sequences2 (db animal-id day)
  (flet ((pos-to-base (pos)
           (consensus-base (apply #'double-to-single-strand
                                  (apply #'filter-bias pos)
                                  ))))
   (loop for (chr-id chr-name) in (sqlite:execute-to-list
                                   db "select id, name from chromosomes;")
      for counts = (sqlite:execute-to-list
                    db "select Af,Ar,Cf,Cr,Gf,Gr,Tf,Tr
                       from pileup where animal = ? and day = ? and chromosome = ?
                       order by position asc;"
                    animal-id day chr-id)
      collecting (make-instance 'seq :name chr-name
                                :characters (map 'vector #'pos-to-base counts)))))

