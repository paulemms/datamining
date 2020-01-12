### Name: mine.associations
### Title: Associations in a contingency table
### Aliases: mine.associations


### ** Examples

data(Titanic)
mine.associations(Titanic)
# Females are twice as likely to survive as the average person.
# Members of 3rd class are twice as likely to be children as the average person.
# Etc.

# focus on associations with survival
mine.associations(Titanic,target="Survived")

# focus on children
mine.associations(Titanic[,,1,],target="Survived")



