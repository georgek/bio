(in-package :gk-bio)

(defun get-consensus-sequences (db animal-id day)
  "Writes consensus from DB for ANIMAL-ID and DAY to STREAM.  DB must be an
  open sqlite database."
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


