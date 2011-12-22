from LinkUtils import LinkCalculator


LINK_FIELDS = ['Source-Prot', 'Target-Prot','Source-Start', 'Source-End', 'Target-Start', 'Target-End',
                'Source-Seq', 'Target-Seq', 'Correct-Num', 'Total-Num', 'This-Score', 'Total-Score', 
                'Source-Cons', 'Target-Cons'] + LinkCalculator().get_fields()
